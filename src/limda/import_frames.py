import pandas as pd
import numpy as np
import pickle
from typing import Union, Any
import pathlib
from typing import Tuple
import yaml
from tqdm import tqdm, trange
from . import const as C
from .import_frame import ImportFrame
from .SimulationFrame import SimulationFrame
import os

class ImportFrames(

):
    """シミュレーションしたデータを読み込むためのクラス
    複数のフレームを読み込む(file -> SimulationFrames)
    """
#-----------------------
    def __init__(self):
        pass
#----------------------------------
    def import_limda_default(self):
        """limdaのデフォルトファイル(.limda.yaml)を読み込む
        """
        limda_dot_path = pathlib.Path.home() / ".limda.yaml"
        if pathlib.Path.exists(limda_dot_path):
            with open(limda_dot_path, "r") as f:
                self.limda_default = yaml.safe_load(f)
#----------------------------------------------------------------------------------------------
    def import_vasp(self, calc_directory: Union[str, pathlib.Path], NELM : int = None):
        """vaspで計算した第一原理MDファイルから、
        原子の座標, cellの大きさ, 原子にかかる力, ポテンシャルエネルギーを読み込む
        Parameters
        ----------
            calc_directory: str
                vaspで計算したディレクトリ
            NELM: int
                最大のIteration回数, 最大のiteration回数に達したframeはimportしない
        Note
        ----
            読み込んだデータ
                simulation_frames[step_idx][['x', 'y', 'z']] : 原子の座標
                simulation_frames[step_idx][['fx', 'fy', 'fz']] : 原子にかかる力
                simulation_frames[step_idx].potential_energy : ポテンシャルエネルギー
                simulation_frames[step_idx].cell : セルサイズ
                simulation_frames[step_idx].virial_tensor : virialテンソル
        """
        if NELM is None:
            if "NELM" in self.limda_default:
                    NELM = self.limda_default["NELM"]
            else:
                NELM = 1e6
        calc_directory = pathlib.Path(calc_directory)
        first_sf = SimulationFrame()
        first_sf.atom_symbol_to_type = self.atom_symbol_to_type
        atom_types = first_sf.import_from_poscar(f'{calc_directory}/POSCAR')

        with open(calc_directory / "OUTCAR", "r") as f:
            lines = f.readlines()
            splines = list(map(lambda l:l.split(), lines))

        for line_idx, spline in enumerate(splines):
            if len(spline) <= 1:
                continue

            if len(spline) == 3 and spline[0] == "POSITION" and spline[1] == "TOTAL-FORCE":
                if iteration == NELM:
                    continue
                sf = SimulationFrame()
                sf.atom_symbol_to_type = self.atom_symbol_to_type
                sf.atom_type_to_mass = self.atom_type_to_mass
                sf.atom_type_to_symbol = self.atom_type_to_symbol 

                sf.cell = np.empty(3, dtype=np.float32)
                for dim in range(3):
                    sf.cell[dim] = float(splines[cell_line_idx+dim+1][dim])

                atoms_dict_keys = ['type','x','y','z','fx','fy','fz']
                atoms_dict = {key: val for key, val in zip(atoms_dict_keys, [[] for i in range(7)])}
                atoms_dict['type'] = atom_types
                for atom_idx in range(len(atom_types)):
                    for key_idx in range(len(atoms_dict_keys)-1):
                        atoms_dict[atoms_dict_keys[key_idx+1]].append(float(splines[line_idx+2+atom_idx][key_idx]))

                sf.atoms = pd.DataFrame(atoms_dict)
                # potential_energy
                sf.potential_energy = float(splines[potential_energy_idx][4])
                # virial tensor
                sf.virial_tensor = np.empty((3, 3), dtype=np.float32)
                for i in range(3):
                    sf.virial_tensor[i][i] = float(splines[virial_tensor_idx][i+1])
                for i in range(3):
                    sf.virial_tensor[i][(i+1)%3] = float(splines[virial_tensor_idx][i+4])
                    sf.virial_tensor[(i+1)%3][i] = float(splines[virial_tensor_idx][i+4])

                self.sf.append(sf)  
 
            if len(splines[line_idx]) == 6 and splines[line_idx][0] == "direct" \
                and splines[line_idx][1] == "lattice":
                cell_line_idx = line_idx
                     
            if len(splines[line_idx]) >= 4 and \
                splines[line_idx][0] == "energy" and \
                splines[line_idx][1] == "without" and \
                splines[line_idx][2] == "entropy":
                potential_energy_idx = line_idx 

            if len(splines[line_idx]) == 7 and splines[line_idx][0] == "Total":
                virial_tensor_idx = line_idx

            if len(splines[line_idx]) == 5 and splines[line_idx][0] == "---------------------------------------" \
                and splines[line_idx][1] == "Iteration":
                iteration = int(splines[line_idx][3][:-1])
#---------------------------------------------------------------------------------------------------
    def import_dumpposes(self, dir_name:Union[str, pathlib.Path]=None, step_nums:list[int]=None, skip_num: int=None):
        """Laichで計算したdumpposを複数読み込む
        Parameters
        ----------
            dir_name: str
                dumpposが入っているフォルダのパス
                指定しないときは、current directryになる
            step_nums: listやイテレータ
                指定したdumpposを読み込む, 
                step_nums=range(0, 301, 100)とすると、
                dump.pos.0, dump.pos.100, dump.pos.200, dump.pos.300を読み込む
            skip_num: int
                いくつおきにdumpposを読み込むのか
                skip_num = 10とすると、10個飛ばしでdumpposを読み込む
        """
        assert self.atom_symbol_to_type is not None, "import atom symbol first"
        assert self.atom_type_to_mass is not None, "import atom symbol first"
        assert self.atom_type_to_symbol is not None, "import atom symbol first" 

        if dir_name is None:
            dir_name = os.getcwd()

        file_names_in_current_dir = os.listdir(dir_name)
        if step_nums is None:
            step_nums = []
            for file_name in file_names_in_current_dir:
                if len(file_name) >= 9 and file_name[:9] == 'dump.pos.':
                    step_nums.append(int(file_name[9:]))

        step_nums.sort()
        if skip_num is not None:
            step_nums = step_nums[::skip_num]

        self.sf = [SimulationFrame() for _ in range(len(step_nums))]
        
        for step_idx, step_num in enumerate(tqdm(step_nums, desc='[importing dumpposes]')):
            self.sf[step_idx].step_num = step_num
            self.sf[step_idx].atom_symbol_to_type = self.atom_symbol_to_type
            self.sf[step_idx].atom_type_to_mass = self.atom_type_to_mass
            self.sf[step_idx].atom_type_to_symbol = self.atom_type_to_symbol
            self.sf[step_idx].import_dumppos(f'{dir_name}/dump.pos.{step_num}')
#-------------------------------------------------------------------------------
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
        if len(atom_symbol_list) == 0:
            assert "para" in self.limda_default
            atom_symbol_list = self.limda_default["para"]
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
#----------------------------------------------------------
    def import_allegro_frames(self, file_name: Union[str, pathlib.Path])->list[dict]:
        """
        pickle fileを読み込み、sfsに入れます。
        pickle file には 
        cell, position, force, atom_types, cut_off, edge_index, potential_energyの情報が入っていますが、
        cut_off, edge_indexの情報はsfs.sfには入らないため、保持するためには返り値であるframesを受け取る必要があります.

        Parameters
        ----------
            file_name :Union[str, pathlib.Path]
                importするpickleファイルのパス
        Return val
        ----------
            frames :list[dict]
                pickelファイルに入っている情報をsfごとにdictで保持したlist
        """
        with open(file_name, 'rb') as p:
            frames = pickle.load(p)
        for frame in frames:
            sf = SimulationFrame()
            sf.cell = frame["cell"]
            sf.potential_energy = frame["potential_energy"]
            sf.atoms = pd.DataFrame(frame["atom_types"] + 1, columns=["type"])
            sf.atoms[["x", "y", "z"]] = pd.DataFrame(frame["pos"])
            sf.atoms[["fx", "fy", "fz"]] = pd.DataFrame(frame["force"])
            sf.virial_tensor = frame["virial"]
            self.sf.append(sf)

        return frames


    