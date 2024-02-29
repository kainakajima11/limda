import pathlib
import os
import numpy as np
import pickle
from tqdm import trange
from .export_frame import ExportFrame
class ExportFrames(
    ExportFrame
):
    def __init__(self):
        pass
#---------------------------------------------------------------------------------
    def export_dumpposes(self, output_folder: str=None, out_columns=None) -> None:
        """SimulationFramesに入ってるSimulationFrameを出力する。
        Parameters
        ----------
        output_folder: str
            出力する場所のパス
        out_columns: list[str]
            dumpposファイルに出力する列の名前の入ったlist
        """
        if output_folder is None:
            output_folder = pathlib.Path.cwd()
        else:
            output_folder = pathlib.Path(output_folder)

        for idx, frame in enumerate(self.sf):
            if frame.step_num is None:
                step_num = idx
            else:
                step_num = frame.step_num
            frame.export_dumppos(ofn = output_folder / f'dump.pos.{step_num}',
                                  time_step = step_num, out_columns=out_columns)
#--------------------------------------------------------------------------------
    def export_allegro_frames(self,
                       output_dir: str,
                       output_file_name: str,
                       cut_off: float,
                       shuffle: bool= False,
                       seed: int=1,
                       test_size: float=None,
                       test_output_dir: str=None,
                       test_output_file_name: str=None,
                       exclude_too_small_cell : bool = True,
                       exclude_too_large_force : bool = True,
                       max_allowable_force : float = 50.0,
                       ):
        """
        allegro用のデータセットを保存する
        Parameters
        ----------
            output_dir : str
                出力する場所
            output_file_name : str
                出力するfile名 {output_file_name}.pickle が出力される
            cut_off : float
                cutoff距離
            shuffle : bool
                フレームをシャッフルするか
            seed : int
                シャッフルするときのシード値
            test_size : float
                test用にする割合
            test_output_dir : str
                test用 : 出力される場所
            test_output_file_name : str
                test用 : 出力されるfile名
            exclude_too_small_cell : bool
                cutoff x 2 以下のセルサイズを持つフレームを除外するか  
            exclude_too_large_force : bool
                forceが基準値(max_allowable_force)より大きいフレームを除外するか
            max_allowable_force : float
                フレームを除外する力の基準値 (exclude_unsuitable_force_frame == True のとき)
        """
        if test_size is not None:
            assert 0.0 <= test_size <= 1.0
            assert test_output_dir is not None
            assert test_output_file_name is not None
            test_output_dir = pathlib.Path(test_output_dir)
            test_output_file_name = pathlib.Path(f"{test_output_file_name}.pickle")
            test_frames_path = test_output_dir / test_output_file_name
            os.makedirs(test_output_dir, exist_ok=True)

        output_dir = pathlib.Path(output_dir)
        output_file_name = pathlib.Path(f"{output_file_name}.pickle")
        frames_path = output_dir / output_file_name
        os.makedirs(output_dir, exist_ok=True)

        train_frames = []
        test_frames = []
        if shuffle:
            self.shuffle_sfs(seed=seed)
        for sf_idx in range(len(self)):
            data = {}
            data["cell"] = np.array(self.sf[sf_idx].cell, dtype=np.float32)
            if  exclude_too_small_cell and np.any(np.diag(self.sf[sf_idx].cell) < 2 * cut_off):
                print(f"Exculuded frame : cellsize(={np.min(np.diag(self.sf[sf_idx].cell))}) is smaller than 2 x cutoff(= {cut_off*2})", flush=True)
                continue
            data["pos"] = np.array(self.sf[sf_idx].atoms[["x","y","z"]].values, dtype=np.float32)
            data["force"] = np.array(self.sf[sf_idx].atoms[["fx","fy","fz"]].values, dtype=np.float32)
            if exclude_too_large_force and np.abs(data['force']).max().item() > max_allowable_force:
                print(f"Exculuded frame : force(={np.abs(data['force']).max().item()}) is larger than reference value of force(={max_allowable_force})", flush=True)
                continue
            data["atom_types"] = np.array(self.sf[sf_idx].atoms["type"].values)
            data["atom_types"] -= 1
            data["cut_off"] = np.array(cut_off, dtype=np.float32)
            data["potential_energy"] = np.array(self.sf[sf_idx].potential_energy, dtype=np.float32)
            data["virial"] = np.array(self.sf[sf_idx].virial_tensor, dtype=np.float32)

            edge_index, shift = self.sf[sf_idx].get_edge_index(cut_off=cut_off)

            data["edge_index"] = np.array(edge_index)
            data["shift"] = np.array(shift, dtype=np.float32)
            if test_size is not None:
                if sf_idx < len(self)*(1.0-test_size):
                    train_frames.append(data)
                else:
                    test_frames.append(data)
            else:
                train_frames.append(data)

        with open(frames_path, "wb") as f:
            pickle.dump(train_frames, f)

        if test_size is not None:
            with open(test_frames_path, "wb") as f:
                pickle.dump(test_frames, f)
#---------------------------------------------------------------------
    def export_lammps_dumpposes(self, ofn: str, out_columns=None) -> None:
            """lammps形式のdumpposを出力する
            Parameters
            ----------
                ofn: str
                    lammps形式のdumpposの出力先
                out_columns: List[str]
                    sdat.atomsのどのカラムを出力するのか
                    デフォルトは['type', 'x', 'y', 'z']
            """
            if out_columns is None:
                out_columns = ['type', 'x', 'y', 'z']

            with open(ofn, 'w') as f:
                f.write('')

            for step_idx in trange(len(self.sf), desc='[exporting lammps dumpposes]'):
                header = []
                header.append(f'ITEM: TIMESTEP\n')
                if self.sf[step_idx].step_num is None:
                    header.append(f'{step_idx}\n')
                else:
                    header.append(f'{self.sf[step_idx].step_num}\n')
                header.append(f'ITEM: NUMBER OF ATOMS\n')
                header.append(f'{self.sf[step_idx].get_total_atoms()}\n')
                header.append(f'ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
                header.append(f'0.0000000000000000e+00 {self.sf[step_idx].cell[0]:.16e} 0.0000000000000000e+00\n')
                header.append(f'0.0000000000000000e+00 {self.sf[step_idx].cell[1]:.16e} 0.0000000000000000e+00\n')
                header.append(f'0.0000000000000000e+00 {self.sf[step_idx].cell[2]:.16e} 0.0000000000000000e+00\n')
                header.append(f'ITEM: ATOMS id {" ".join(out_columns)}\n')

                with open(ofn, 'a') as f:
                    f.writelines(header)
                
                # 1-index
                self.sf[step_idx].atoms.index += 1
                self.sf[step_idx].atoms.to_csv(ofn, columns=out_columns, sep=' ', header=None, mode='a')
                # 0-index
                self.sf[step_idx].atoms.index -= 1

    #--------------------------------------------------------------------------------
    def export_xyzs(self, output_folder: str, out_columns=None) -> None:
        output_folder = pathlib.Path(output_folder)
        output_folder.mkdir(exist_ok=True, parents=True)
        for frame_idx in range(len(self.sf)):
            self.sf[frame_idx].export_xyz(output_folder / f"{frame_idx}.xyz", out_columns=out_columns)


