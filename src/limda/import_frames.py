import pandas as pd
import numpy as np
from typing import Union
import pathlib
from typing import Tuple
from tqdm import tqdm, trange
import limda.const as C
from .import_frame import ImportFrame
from .SimulationFrame import SimulationFrame
import os

class ImportFrames(
    ImportFrame
):
    """シミュレーションしたデータを読み込むためのクラス
    複数のフレームを読み込む(file -> SimulationFrames)
    """
#-----------------------
    def __init__(self):
        pass
#----------------------------------------------------------------------------------------------
    def import_vasp(self, calc_directory: Union[str, pathlib.Path]):
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
        calc_directory = pathlib.Path(calc_directory)
        first_sf = SimulationFrame()
        first_sf.atom_symbol_to_type = self.atom_symbol_to_type
        atom_types = first_sf.import_from_poscar(f'{calc_directory}/POSCAR')

        with open(calc_directory / "OUTCAR", "r") as f:
            lines = f.readlines()
            splines = list(map(lambda l:l.split(), lines))

        for line_idx, spline in enumerate(splines):
            if len(spline) == 0:
                continue
            if len(spline) == 3 and spline[0] == "POSITION" and spline[1] == "TOTAL-FORCE":
                sf = SimulationFrame()
                sf.atom_symbol_to_type = self.atom_symbol_to_type
                sf.atom_type_to_mass = self.atom_type_to_mass
                sf.atom_type_to_symbol = self.atom_type_to_symbol 

                sf.cell = [None, None, None]
                for dim in range(3):
                    sf.cell[dim] = float(splines[cell_line_idx+dim+1][dim])

                atoms_dict_keys = ['type','x','y','z','fx','fy','fz']
                atoms_dict = {key: val for key, val in zip(atoms_dict_keys, [[] for i in range(7)])}
                atoms_dict['type'] = atom_types
                for atom_idx in range(len(atom_types)):
                    for key_idx in range(len(atoms_dict_keys)-1):
                        atoms_dict[atoms_dict_keys[key_idx+1]].append(float(splines[line_idx+2+atom_idx][key_idx]))

                sf.atoms = pd.DataFrame(atoms_dict)
                sf.potential_energy = float(splines[potential_energy_idx][4])
                self.sf.append(sf)  
 
            if len(splines[line_idx]) == 6 and splines[line_idx][0] == "direct" \
                and splines[line_idx][1] == "lattice":
                cell_line_idx = line_idx
                     
            if len(splines[line_idx]) >= 4 and \
                splines[line_idx][0] == "energy" and \
                splines[line_idx][1] == "without" and \
                splines[line_idx][2] == "entropy":
                potential_energy_idx = line_idx 
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

        for step_idx, step_num in enumerate(tqdm(step_nums)):
            self.sf[step_idx].step_num = step_num
            self.sf[step_idx].atom_symbol_to_type = self.atom_symbol_to_type
            self.sf[step_idx].atom_type_to_mass = self.atom_type_to_mass
            self.sf[step_idx].atom_type_to_symbol = self.atom_type_to_symbol
            self.sf[step_idx].import_dumppos(f'{dir_name}/dump.pos.{step_num}')
                        
                
     


    