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
#-------------------------------------------------------------------------------
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