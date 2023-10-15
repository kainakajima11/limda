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
                       ):
        """allegro用のデータセットを保存する
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
            data["pos"] = np.array(self.sf[sf_idx].atoms[["x","y","z"]].values, dtype=np.float32)
            data["force"] = np.array(self.sf[sf_idx].atoms[["fx","fy","fz"]].values, dtype=np.float32)
            data["atom_types"] = np.array(self.sf[sf_idx].atoms["type"].values)
            data["atom_types"] -= 1
            data["cut_off"] = np.array(cut_off, dtype=np.float32)
            data["potential_energy"] = np.array(self.sf[sf_idx].potential_energy, dtype=np.float32)
            data["virial"] = np.array(self.sf[sf_idx].virial_tensor, dtype=np.float32)

            edge_index = [[],[]]
            edge_index = self.sf[sf_idx].get_edge_idx(cut_off=cut_off)

            data["edge_index"] = np.array(edge_index)
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
    def concat_simulation_frames(self, sfs_list:list):
        """sfsを結合する
        Parameters
        ----------
            sfs_list : list[SimulationFrames]
                結合するsfsのリスト, 
        Note
        ----
            concat_sfsメソッドを使用するSimulationFramesは
            import_para()後のを使う
        """
        self.sf = []
        for outer_sfs in simulation_frames_list:
            self.sf.extend(outer_sfs.sf)
        
        self.step_nums = list(range(len(self.sdat)))
        self.step_num_to_step_idx = {
            step_num:step_idx for step_idx, step_num in enumerate(range(len(self.step_nums)))
        }
