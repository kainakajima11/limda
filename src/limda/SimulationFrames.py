import pandas as pd
import numpy as np
import pathlib
import random
import os
from tqdm import tqdm, trange
from .import_frames import ImportFrames
from .export_frames import ExportFrames
from .SimulationFrame import SimulationFrame

class SimulationFrames(
    ImportFrames,
    ExportFrames
):
    """シミュレーションしたデータを読み込み、書き込み、分析するためのクラス
    複数のフレームを同時に扱う

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
    def shuffle_sfs(self, seed:int=1):
        """self.sfの順番をシャッフルする
        Parameters
        ----------
            seed: int
                乱数seed値
        """
        random.seed(seed)
        random.shuffle(self.sf)
# -----------------------------------------------------------------  
    def concat_sfs(self, simulation_frames_list:list):
        """sfsを結合する
        Parameters
        ----------
            simulation_frames_list : list[SimulationFrames]
                結合するsfsのリスト, 
        Note
        ----
            concat_sfsメソッドを使用するSimulationFramesは
            import_para()後のを使う
        """
        self.sf = []
        for outer_sfs in simulation_frames_list:
            for step_idx in range(len(outer_sfs)):
                self.sf.append(outer_sfs.sf[step_idx])
        
        step_nums = list(range(len(self.sf)))
        for step_idx, step_num in enumerate(tqdm(step_nums)):
            self.sf[step_idx].step_num = step_num
#-------------------------------------------------------------------
    def split_sfs_specified_list_size(self, list_size: int)->list:
        """sfsを複数のsfsに分け, sfsのlistを返す。
            listのサイズを指定できる
        Parameters
        ----------
            list_size: int
                sfsを何分割するか
        Return val
        ----------
            sfs_list: list[SimulationFrames()]
            元のsfsを複数に分けたときのsfsから成るlist
        Example
        -------
            sfs  = SimulationFrames()  # len(sfs) = 10
            sfs_list = sfs.split_sframes_specify_list_size(3)
                ->sfs_listは3つのsfsからなるlistで、
                    len(sfs_list[i]) = [4,3,3]
        """
        sfs_list = [SimulationFrames() for _ in range(list_size)]
        item_num_list = [int(len(self)/list_size) for _ in range(list_size)]
        for i in range(len(self) % list_size):
            item_num_list[i] += 1
        for idx, item in enumerate(item_num_list):
            for i in range(1, item+1):
                sfs_list[idx].sf.append(i+sum(item_num_list[0:idx]))
        for frames in sfs_list:
            frames.atom_symbol_to_type = self.atom_symbol_to_type
            frames.atom_type_to_symbol = self.atom_type_to_symbol
            frames.atom_type_to_mass = self.atom_type_to_mass
        return sfs_list
#---------------------------------------------------------------------------------------------------------------
    def split_sfs(self, each_sfs_size: int, keep_remains:bool = False)->list:
        """ sfsを複数のsfsに分け, sfsのlistを返す。
            sfs1つ1つサイズを指定できる。
        Parameters
        ----------
            each_sfs_size: int
                list内のsfsのサイズ
            keep_remains: bool
                残りを捨てるか、後ろにくっつけるか
        Return val
        ----------
            sfs_list: list[SimulationFrames()]
            元のsfsを複数に分けたときのsfsから成るlist
        Example
        -------
            sfs  = SimulationFrames()  # len(sfs) = 10
            sfs_list = sfs.split_sframes(3)
                -> len(sfs_list[i]) = [3,3,3] (keep_remains = False)
                len(sfs_list[i]) = [3,3,3,1] (keep_remains = True)
        """
        list_size = int(len(self) / each_sfs_size)
        main_sfs = SimulationFrames()
        remain_sfs = SimulationFrames()
        remain = len(self)%each_sfs_size
        main_sfs.sf = self.sf[0:len(self)-remain]
        remain_sfs.sf = self.sf[len(self)-remain:len(self)]
        sfs_list = main_sfs.split_sframes_specify_list_size(list_size)
        if keep_remains and remain != 0:
            remain_sfs.atom_symbol_to_type = self.atom_symbol_to_type
            remain_sfs.atom_type_to_symbol = self.atom_type_to_symbol
            remain_sfs.atom_type_to_mass = self.atom_type_to_mass
            sfs_list.append(remain_sfs)
        return sfs_list
