import pathlib
import os
import numpy as np
import pickle
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
            self.shuffle_sf(seed=seed)
        for sf_idx in range(len(self)):
            data = {}
            data["cell"] = np.array(self.sf[sf_idx].cell, dtype=np.float32)
            data["pos"] = np.array(self.sf[sf_idx].atoms[["x","y","z"]].values, dtype=np.float32)
            data["force"] = np.array(self.sf[sf_idx].atoms[["fx","fy","fz"]].values, dtype=np.float32)
            data["atom_types"] = np.array(self.sf[sf_idx].atoms["type"].values)
            data["atom_types"] -= 1
            data["cut_off"] = np.array(cut_off, dtype=np.float32)
            data["potential_energy"] = np.array(self.sf[sf_idx].potential_energy, dtype=np.float32)

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
