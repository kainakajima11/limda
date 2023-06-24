import pandas as pd
import numpy as np
from typing import Union
import pathlib 
import re
import sys
import limda.const as C

class ImportFrame(

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
#----------------------
    def __init__(self):
        pass
#-----------------------------------------------------------------------   
    def import_input(self, file_path: Union[str, pathlib.Path]) -> None: #ky
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
#---------------------------------------------------------------
    def import_para_from_list(self, atom_symbol_list:list[str]):
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
        for atom_symbol, atom_type in self.atom_symbol_to_type.items():
            self.atom_type_to_mass[atom_type] = C.ATOM_SYMBOL_TO_MASS[atom_symbol]
#-------------------------------------------------------
    def import_para_from_str(self, atom_symbol_str:str):
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
        self.import_para_from_list(atom_symbol_str.split())
#---------------------------------------------------------------------
    def import_car(self, file_path: Union[str, pathlib.Path]) -> None:
        """ Car file を読み込み 
            self.atoms, self.cell の更新
            Parameter
            ----------
                file_path: Union[str, Path]
                carfileのpath
        """
        car_df = pd.read_csv(file_path, names=['symbol+id', 'x', 'y', 'z', 'XXXX', '1', 'xx', 'symbol', '0.000'],
                         usecols=['x', 'y', 'z', 'symbol'],
                         skiprows=4,  sep='\s+')
        self.cell = np.float_(car_df.iloc[0, 0:3]) # cell size部分を抜き取る
        car_df = car_df[1:].dropna() 
        car_df.insert(0, 'type', car_df['symbol'].map(self.atom_symbol_to_type)) # type列を作成
        self.atoms = car_df[['type', 'x', 'y', 'z']] 
#-------------------------------------------------------------------------
    def import_dumppos(self, file_path: Union[str, pathlib.Path]) -> None: 
        """ dumppos file を読み込み 
            self.atomsの更新
             Parameter
            ----------
                file_path: Union[str, Path]
                dumpposfileのpath
        """
        current_row = 0
        with open(file_path, 'r') as ifp:
            while True:
                current_row += 1
                spline = ifp.readline().split()
                if len(spline) == 0:
                    continue
                if spline[0] == "ITEM:" and spline[1] == "BOX":
                    self.cell = [None, None, None]
                    for dim in range(3):
                        spline = ifp.readline().split()
                        self.cell[dim] = float(spline[1])
                    current_row += 3
                    continue
                if spline[0] == "ITEM:" and spline[1] == 'ATOMS':
                    columns = spline[3:]
                    break

        self.atoms = pd.read_csv(
            file_path, skiprows = current_row, sep='\s+', names=columns)
        if 'type' in self.atoms:
            self.atoms['type'] = self.atoms['type'].astype(int)
        if 'mask' in self.atoms:
            self.atoms['mask'] = self.atoms['mask'].astype(int)

        self.atoms.index = self.atoms.index - 1
        self.atoms.sort_index(inplace=True)
#------------------------------------------
    def import_mol(self, molecular_fomula): #ky
            '''分子式から原子配置を読み込む
            Parameters
            ----------
                molecular_fomula : str
                    読み込む分子式
            '''
            try:
                from ase import build
            except:
                pass

            ase_atoms = build.molecule(molecular_fomula)
            self.atoms = pd.DataFrame(ase_atoms.positions, columns=['x', 'y', 'z'])
            self.atoms['type'] = ase_atoms.get_chemical_symbols()
            self.atoms['type'] = self.atoms['type'].map(self.atom_symbol_to_type)
            self.atoms['x'] += abs(self.atoms['x'].min()) + 0.1 #?
            self.atoms['y'] += abs(self.atoms['y'].min()) + 0.1
            self.atoms['z'] += abs(self.atoms['z'].min()) + 0.1
#------------------------------------------------------------------------------------
    def import_from_poscar(self, poscar_path: Union[str, pathlib.Path]) -> list[int]:
        """vaspに用いるPOSCARから, 
        原子それぞれの種類を表すリストを作成する。
        また、初期構造のSimulationFrame(原子の座標のみ)が得られる。
        Parameters
        ----------
            poscar_path: Union[str, Path]
                vaspで計算したディレクトリ内のPOSCARのpath 
        Return val
        ----------
            atom_types: list[int]
            原子の種類をtype listと照らし合した時の整数が入っています。  

            sf: SimulationFrame
            t=0 の SimulationFrame
        """
        with open(poscar_path, "r") as f:
            f.readlines(2)
            self.cell = np.array([None, None, None])
            for dim in range(3):
                self.cell[dim] = float(f.readline().split()[dim])
            atom_symbol_list = list(f.readline().split())
            atom_type_counter = list(map(int, f.readline().split()))
            atom_types = []
            for atom_type_count, atom_symbol in zip(atom_type_counter, atom_symbol_list):
                for _ in range(atom_type_count):
                    atom_types.append(self.atom_symbol_to_type[atom_symbol])

            self.atoms = pd.read_csv(
                f, skiprows = 1, sep='\s+', names=("x", "y", "z"))
            
        return atom_types
#-----------------------------------------------------------------------------------
    def import_xyz(self, ifn: Union[str,pathlib.Path])->None:
        """xyz fileを読み込む.
        Parameter
        ---------
        ifn: Union[str,Path]
            読み込むfile名
        """
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()
        total_atom = int(lines[0])
        lattice_value = re.search('Lattice="(.*?)"', lines[1])

        if lattice_value is not None:
            cellsize = lattice_value.group(1).split()
            self.cell[0] = float(cellsize[0])
            self.cell[1] = float(cellsize[1])
            self.cell[2] = float(cellsize[2])
        
        splines = np.array([l.split() for l in lines[2:2+total_atom]])
        atom_data = dict()
        try: 
            atom_data['type'] = splines[:, 0].astype(int)
        except:
            atom_symbols = splines[:,0].astype(str)
            if self.atom_symbol_to_type is None:
                print("error : atom_symbol_to_type is not defined")
                print("Import para first")
                sys.exit(-1)
            atom_data['type'] = np.array(
                [self.atom_symbol_to_type[atom_symbol] for atom_symbol in atom_symbols]
            )

        for idx, dim in enumerate(['x', 'y', 'z']):
            atom_data[dim] = splines[:, idx+1].astype(float)

        index = np.arange(total_atom)
        self.atoms = pd.DataFrame(data=atom_data, index=index)
        self.atoms.sort_index(inplace=True)