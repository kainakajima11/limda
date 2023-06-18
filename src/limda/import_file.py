import pandas as pd
import numpy as np
from typing import Union
from pathlib import Pathgit 
from limda.import_atoms_symbol_to_mass import import_atoms_symbol_to_mass
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

    def import_atom_from_list(self, atom_symbol_list):
        """原子のリストからatom_symbol_to_type, atom_type_to_symbol, atom_type_to_massを作成する.
        Parameters
        ----------
            atom_symbol_list : list
                原子のリスト

        Example
        -------
            atom_symbol_list = ['C', 'H', 'O', 'N']
            の場合、Cの原子のタイプが1, Hの原子のタイプが2, Oの原子のタイプが3, Nの原子のタイプが4となる

        """ 
        atom_symbol_to_type = {}
        type_list = [i for i in range(1, len(atom_symbol_list)+1)]
        atom_symbol_to_type = {key: val for key, val in zip(atom_symbol_list, type_list)}
        # type -> symbol# symbol -> type # type -> mass
        self.atom_symbol_to_type = atom_symbol_to_type
        self.atom_type_to_symbol = {
            atom_type: atom_symbol for atom_symbol, atom_type in self.atom_symbol_to_type.items()}
        self.atom_type_to_mass = {}
        atom_symbol_to_mass = import_atoms_symbol_to_mass()
        for atom_symbol, atom_type in self.atom_symbol_to_type.items():
            self.atom_type_to_mass[atom_type] = atom_symbol_to_mass[atom_symbol]

    
    def import_atom_from_str(self, atom_symbol_str):
        """
            受け取ったstrをlistにして、
            import_para_from_list()を呼び出す。
            Parameters
            ----------
                para_atom_symbol_list : list   
                空白区切りの原子の文字列

            Example
            -------
                para_atom_symbol_str = 'C H O N' #原子と原子の間には、スペース
                の場合、Cの原子のタイプが1, Hの原子のタイプが2, Oの原子のタイプが3, Nの原子のタイプが4となる
        """
        self.import_atom_from_list(atom_symbol_str.split())



    def import_car(self, file_path: Union[str, Path]) -> None:
        """ Car file を読み込み 
            self.atoms, self.cell の更新
        """
        car_df = pd.read_csv(file_path, names=['symbol+id', 'x', 'y', 'z', 'XXXX', '1', 'xx', 'symbol', '0.000'],
                         usecols=['x', 'y', 'z', 'symbol'],
                         skiprows=4,  sep='\s+')
        self.cell = np.float_(car_df.iloc[0, 0:3]) # cell size部分を抜き取る
        car_df = car_df[1:].dropna() 
        car_df.insert(0, 'type', car_df['symbol'].map(self.atom_symbol_to_type)) # type列を作成
        self.atoms = car_df[['type', 'x', 'y', 'z']] 
        