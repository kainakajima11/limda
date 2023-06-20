import pandas as pd
import numpy as np
import pathlib
import random
import os
from tqdm import tqdm, trange
from .import_frame import ImportFrame
from .import_frames import ImportFrames
from .SimulationFrame import SimulationFrame

class SimulationFrames(
    ImportFrames,
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
        self.sf:list[SimulationFrame] = []
        self.atom_symbol_to_type: dict[str, int] = None
        self.atom_type_to_symbol : dict[int, str] = None
        self.atom_type_to_mass : dict[int, float] = None
        pass
#---------------------
    def __len__(self):
        """
        len(sfs)でlen(sfs.sf)を得ることができる
        """
        return len(self.sf)
#----------------------------- 
    def __getitem__(self, key):
        """sfs = SimulationFrames()
        sfs[step_idx]でsfs.sdat[step_idx]を得ることができる
        """
        return self.sf[key]
#------------------------------------
    def shuffle_sf(self, seed:int=1):
        """self.sfの順番をシャッフルする
        Parameters
        ----------
            seed: int
                乱数seed値
        """
        random.seed(seed)
        random.shuffle(self.sf)
# -----------------------------------------------------------------  
    def concat_simulation_frames(self, simulation_frames_list:list):
        """sfsを結合する
        Parameters
        ----------
            simulation_frames_list : List[SimulationFrames]
                結合するsfsのリスト, 
        Note
        ----
            concat_sdatsメソッドを使用するSimulationFramesは
            import_para()後のを使う
        """
        self.sf = []
        for outer_sfs in simulation_frames_list:
            for step_idx in range(len(outer_sfs)):
                self.sf.append(outer_sfs.sf[step_idx])
        
        step_nums = list(range(len(self.sf)))
        for step_idx, step_num in enumerate(tqdm(step_nums)):
            self.sf[step_idx].step_num = step_num

