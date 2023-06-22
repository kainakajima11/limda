import pandas as pd
import numpy as np
import sys
import os
from .import_frame import ImportFrame

class SimulationFrame(
    ImportFrame,
):
    """シミュレーションしたデータを読み込み、書き込み、分析するためのクラス
    一つのフレームを扱う

    Attributes
    ----------
    atoms : pd.DataFrame
        原子のtype, 座標, 速度, 加速度, 力, type, 電荷などを含むpandasのDataFrame
        原子のタイプは1-indexed
    cell : np.array
        cellの大きさが入ったarray, shape:[3]
        cell[0]:x方向, cell[1]:y方向, cell[2]:z方向
    atom_symbol_to_type : dict[str, int]
        原子のシンボルをkey, 原子のtypeをvalueとするdict
    atom_type_to_symbol : dict[int, str]
        原子のtypeをkey, 原子のシンボルをvalueとするdict
    atom_type_to_mass : dict[int, float]
        原子のtypeをkey, 原子の質量(g/mol)をvalueとするdict
    step_num: int
        MDした時のこのSimulationFrameが持つステップ数
    potential_energy: int
        このフレームが持つポテンシャルエネルギー, 単位はeV
    """
    atoms: pd.DataFrame
    cell: np.array # shape:[3]
    atom_symbol_to_type: dict[str, int]
    atom_type_to_symbol : dict[int, str]
    atom_type_to_mass : dict[int, float]
    step_num: int
#----------------------
    def __init__(self):
        self.atoms = None
        self.cell = None
        self.atom_symbol_to_type = None
        self.atom_type_to_symbol = None
        self.atom_type_to_mass = None
        self.step_num = None
        self.potential_energy = None
#-------------------------------------
    def __getitem__(self, key) -> pd.DataFrame:
        """
        sdat.atoms[column]をsdat[column]と省略して書くことが出来る。
        """
        return self.atoms[key]
#-------------------------------------------
    def __setitem__(self, key, val) -> None:
        self.atoms[key] = val
#------------------------------
    def __len__(self) -> int:
            """
            sdat.get_total_atoms()をlen(sdat)と省略して書くことが出来る。
            """
            return self.get_total_atoms()

#-------------------------------------
    def get_total_atoms(self) -> int:
        """
        全原子数を返す関数。
        """
        assert self.atoms is not None, 'Import file first'
        return len(self.atoms)

#--------------------------------------
    def get_atom_type_set(self) -> set:
        """
        系内の原子のtypeのsetを返す関数
        """
        return set(self.atoms['type'])
#----------------------------------------------
    def wrap_atoms(self) -> None: #ky
        """
        セルの外にはみ出している原子をセルの中に入れる。
        """
        assert self.cell is not None, "set sf.cell"
        assert 0 not in set(self.cell), "cell size must not be 0" 

        self.atoms[['x', 'y', 'z']] %= self.cell
#---------------------------------------------------------------------------------
    def replicate_atoms(self, replicate_directions:list[int] = [1, 1, 1]) -> None:
        """
        x, y, z 方向にセルを複製する関数

        Parameters
        ----------
        replicate_directions : list
            x, y, z方向に何倍するかを指定する。
            例えばx方向に2倍,y方向に3倍,z方向に4倍したい時は
            replicate_directions = [2, 3, 4] とする

        """
        for i, dim in enumerate(['x', 'y', 'z']):
            copy_atoms = self.atoms.copy()
            for idx in range(1,replicate_directions[i]): 
                replicate = lambda d: d[dim] + self.cell[i]*idx
                append_atoms = copy_atoms.copy()
                append_atoms[dim] = append_atoms.apply(replicate,axis=1)
                self.atoms = pd.concat([self.atoms, append_atoms])

        self.atoms.reset_index(drop=True, inplace=True)
        for dim in range(3):
            self.cell[dim] *= replicate_directions[dim]
#------------------------------------------------------
    def concat_atoms(self, outer_sf) -> None:
        """
        sfとouter_sfを結合する関数

        Parameters
        ----------
        outer_sf : SimulationFrame
            取り入れたいsfを指定する。

        """
        self.atoms = pd.concat([self.atoms, outer_sf.atoms])
        self.atoms.reset_index(drop=True, inplace=True)
        for dim in range(3):
            self.cell[dim] = max(self.cell[dim], outer_sf.cell[dim])

#------------------------------------------------
    def delete_atoms(self, condition, reindex):
            """
            条件に当てはまる原子を削除する

            Parameters
            ----------
            condition : function
                削除したい原子の条件を指定する関数
            reindex : bool
                reindex == Trueの時は原子のid(idx)は新しく割り振られる
                reindex == Falseの時は原子のid(idx)は新しく割り振られず、削除前のものと変わらない
            """
            if callable(condition):
                target_atoms = condition(self)
                self.atoms = self.atoms[~target_atoms]
            else:
                self.atoms = self.atoms[~condition]
            if reindex:
                self.atoms.reset_index(drop=True, inplace=True)



