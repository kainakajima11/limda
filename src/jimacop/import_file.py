import pandas as pd
import numpy as np
from typing import Union
from pathlib import Path

class ImportFile(

):
    """シミュレーションしたデータを読み込むためのクラス
    一つのフレームを読み込む
    """
    atoms: pd.DataFrame
    cell: np.array # shape:[3]
    atom_symbol_to_type: dict[str, int]
    atom_type_to_symbol : dict[int, str]
    atom_type_to_mass : dict[int, float]
    step_num: int

    def __init__(self):
        pass

    
    def import_input(self, file_path: Union[str, Path]) -> None:
        """Laichのinputファイルを読み込み、
        sf.cell, sf.atomsを更新する
        Parameters
        ----------
            file_path: Union[str, Path]
                input.rdのパス
        Note
        ----
            #fix, #move, #press, #strain, #connect, #thermofree
            #wall, #molecule,は読み込まれません
        """
        self.cell = np.array([0.0, 0.0, 0.0])
        
        if self.atom_symbol_to_type is None:
            raise RuntimeError("Import para first")

        with open(file_path, 'r') as ifp:
            lines = ifp.readlines()

        for idx, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "#cellx":
                self.cell[0] = float(spline[2])

            if spline[0] == "#celly":
                self.cell[1] = float(spline[2])

            if spline[0] == "#cellz":
                self.cell[2] = float(spline[2])

            if spline[0] == "#masses":
                elem_num = int(spline[1])
                for _line in lines[idx+1:idx+1+elem_num]:
                    _spline = _line.split()
                    self.atom_type_to_mass[int(_spline[0])] = float(_spline[1])

            atom_data = dict()
            if spline[0] == "#atoms":
                splines = np.array([l.split()
                                   for l in lines[idx+1:idx+1+int(spline[1])]])
                # 0-indexed
                index = splines[:, 0].astype(int) - 1
                atom_data['type'] = splines[:, 1].astype(int)
                atom_data['mask'] = splines[:, 2].astype(int)
                atom_data['x'] = splines[:, 3].astype(float)
                atom_data['y'] = splines[:, 4].astype(float)
                atom_data['z'] = splines[:, 5].astype(float)
                try:
                    atom_data['vx'] = splines[:, 6].astype(float)
                    atom_data['vy'] = splines[:, 7].astype(float)
                    atom_data['vz'] = splines[:, 8].astype(float)
                except:
                    pass

                self.atoms = pd.DataFrame(data=atom_data, index=index)

