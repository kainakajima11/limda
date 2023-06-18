import pandas as pd
import numpy as np
import pathlib
import random
import os
from tqdm import tqdm, trange
from .import_file import ImportFile
from .SimulationFrame import SimulationFrame

class SimulationFrames(
    ImportFile
):
    """シミュレーションしたデータを読み込み、書き込み、分析するためのクラス
    複数ののフレームを同時に扱う

    Attributes
    ----------
    sf : list[SimulationFrame]
        シミュレーションしたデータを読み込み、書き込み、分析するためのクラス
    atom_symbol_to_type : dict[str, int]
        原子のシンボルをkey, 原子のtypeをvalueとするdict
    atom_type_to_symbol : dict[int, str]
        原子のtypeをkey, 原子のシンボルをvalueとするdict
    atom_type_to_mass : dict[int, float]
        原子のtypeをkey, 原子の質量(g/mol)をvalueとするdict
    
    """
    sf: list[SimulationFrame]
    atom_symbol_to_type: dict[str, int]
    atom_type_to_symbol : dict[int, str]
    atom_type_to_mass : dict[int, float]
#----------------------
    def __init__(self):
        pass
#---------------------
    def __len__(self):
        return len(self.step_nums)
#-----------------------------
    def __getitem__(self, key):
        """sfs = SimulationFrames()
        sfs[step_idx]でsfs.sdat[step_idx]を得ることができる
        """
        return self.sdat[key]
#----------------------------------------------------------    
    def import_dumpposes(self, dir_name:str=None, step_nums:list[int]=None, skip_num: int=None):
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
            self.step_nums = []
            for file_name in file_names_in_current_dir:
                if len(file_name) >= 9 and file_name[:9] == 'dump.pos.':
                    self.step_nums.append(int(file_name[9:]))

        else:
            self.step_nums = []
            for step_num in step_nums:
                self.step_nums.append(step_num)
        self.step_nums.sort()
        if skip_num is not None:
            self.step_nums = self.step_nums[::skip_num]

        self.step_num_to_step_idx = {
            step_num: step_idx for step_idx, step_num in enumerate(self.step_nums)
        }
        self.sdat = [SimulationFrame() for _ in range(len(self.step_nums))]

        for step_idx, step_num in enumerate(tqdm(self.step_nums)):
            self.sdat[step_idx].atom_symbol_to_type = self.atom_symbol_to_type
            self.sdat[step_idx].atom_type_to_mass = self.atom_type_to_mass
            self.sdat[step_idx].atom_type_to_symbol = self.atom_type_to_symbol
            self.sdat[step_idx].import_dumppos(f'{dir_name}/dump.pos.{step_num}')
#--------------------------------------------------------------------------------
    def import_vasp(self, calc_directory: str):
        """vaspで計算した第一原理MDファイルから、
        原子の座標, cellの大きさ, 原子にかかる力, ポテンシャルエネルギーを読み込む
        Parameters
        ----------
            calc_directory: str
                vaspで計算したディレクトリ
        Note
        ----
            読み込んだデータ
                simulation_frames[step_idx][['x', 'y', 'z']] : 原子の座標
                simulation_frames[step_idx][['fx', 'fy', 'fz']] : 原子にかかる力
                simulation_frames[step_idx].potential_energy : ポテンシャルエネルギー
                simulation_frames[step_idx].cell : セルの大きさ
        """
        self.sdat = []
        calc_directory = pathlib.Path(calc_directory)
        with open(calc_directory / "POSCAR", "r") as f:
            for _ in range(5):
                f.readline()
            atom_symbol_list = list(f.readline().split())
            atom_type_counter = list(map(int, f.readline().split()))
            atom_types = []
            for atom_type_count, atom_symbol in zip(atom_type_counter, atom_symbol_list):
                for _ in range(atom_type_count):
                    atom_types.append(self.atom_symbol_to_type[atom_symbol])

        with open(calc_directory / "OUTCAR", "r") as f:
            lines = f.readlines()
            splines = list(map(lambda l:l.split(), lines))

        
        for line_idx, spline in enumerate(splines):
            if len(spline) == 0:
                continue
            if len(spline) == 3 and spline[0] == "POSITION" and spline[1] == "TOTAL-FORCE":
                position_line_idx = line_idx
                sf = SimulationFrame()
                sf.atom_symbol_to_type = self.atom_symbol_to_type
                sf.atom_type_to_mass = self.atom_type_to_mass
                sf.atom_type_to_symbol = self.atom_type_to_symbol 
                cell_line_idx = position_line_idx
                while True:
                    if len(splines[cell_line_idx]) == 6 and splines[cell_line_idx][0] == "direct" \
                        and splines[cell_line_idx][1] == "lattice":
                        break
                    cell_line_idx -= 1
                potential_energy_idx = position_line_idx
                while True:
                    if len(splines[potential_energy_idx]) >= 4 and \
                        splines[potential_energy_idx][0] == "energy" and \
                        splines[potential_energy_idx][1] == "without" and \
                        splines[potential_energy_idx][2] == "entropy":
                        break
                    potential_energy_idx -= 1
                        
                
                sf.cell = [None, None, None]
                sf.cell[0] = float(splines[cell_line_idx+1][0])
                sf.cell[1] = float(splines[cell_line_idx+1+1][1])
                sf.cell[2] = float(splines[cell_line_idx+1+2][2])
                atoms_dict = {
                    'type':atom_types,
                    'x':[],
                    'y':[],
                    'z':[],
                    'fx':[],
                    'fy':[],
                    'fz':[],
                              }
                for atom_idx in range(len(atom_types)):
                    atoms_dict['x'].append(float(splines[position_line_idx+2+atom_idx][0]))
                    atoms_dict['y'].append(float(splines[position_line_idx+2+atom_idx][1]))
                    atoms_dict['z'].append(float(splines[position_line_idx+2+atom_idx][2]))
                    atoms_dict['fx'].append(float(splines[position_line_idx+2+atom_idx][3]))
                    atoms_dict['fy'].append(float(splines[position_line_idx+2+atom_idx][4]))
                    atoms_dict['fz'].append(float(splines[position_line_idx+2+atom_idx][5]))
                sf.atoms = pd.DataFrame(atoms_dict)
                sf.potential_energy = float(splines[potential_energy_idx][4])
                self.sdat.append(sf)
        self.step_nums = list(range(1, len(self.sf) + 1))
        self.step_num_to_step_idx = {
            step_num: step_idx for step_idx, step_num in enumerate(self.step_nums)
        }
#--------------------------------------
    def shuffle_sf(self, seed:int=1):
        """SimulationFrames.sfの順番をシャッフルする
        Parameters
        ----------
            seed: int
                乱数seed値
        """
        random.seed(seed)
        random.shuffle(self.sf)



