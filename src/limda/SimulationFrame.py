import pandas as pd
import numpy as np
import random
import sys
import os
from .import_frame import ImportFrame
from .export_frame import ExportFrame
from .calculate import Calculate
import limda.const as C

class SimulationFrame(
    ImportFrame,
    ExportFrame,
    Calculate
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
#-------------------------------------------------------------------------------------------
    def density(self, x_max=None, x_min=None, y_max=None, y_min=None,z_max=None,z_min=None): #k
        """セル内の密度を計算する関数
        Parameters
        ----------
        x_max: float
            xの上限
        x_min: float
            xの下限
        y_max: float
            yの上限
        y_min: float
            yの下限
        z_max: float
            zの上限
        z_min: float
            zの下限
        Returns
        -------
        density : float
            セル内の密度(g/cm^3)
        """
        x_mx = x_max if x_max is not None else self.cell[0]
        x_mn = x_min if x_min is not None else 0
        assert x_mx >= x_mn, 'set correct x_max and x_min'

        y_mx = y_max if y_max is not None else self.cell[1]
        y_mn = y_min if y_min is not None else 0
        assert y_mx >= y_mn, 'set correct y_max and y_min'

        z_mx = z_max if z_max is not None else self.cell[2]
        z_mn = z_min if z_min is not None else 0
        assert z_mx >= z_mn, 'set correct z_max and z_min'

        # 体積(cm^3)
        volume = (x_mx - x_mn) * (y_mx - y_mn) * (z_mx - z_mn) * (10 ** - 24)
        all_weight = 0
        def condition(sf):
            target_atoms = (x_mn <= sf.atoms['x'])&(sf.atoms['x'] <= x_mx)
            target_atoms &= (y_mn <= sf.atoms['y'])&(sf.atoms['y'] <= y_mx)
            target_atoms &= (z_mn <= sf.atoms['z'])&(sf.atoms['z'] <= z_mx)
            return target_atoms
        atom_type_counter = self.count_atom_types(res_type='dict', condition=condition)
        for atom_symbol, atom_type_count in atom_type_counter.items():
            atom_type = self.atom_symbol_to_type[atom_symbol]
            atom_mass = self.atom_type_to_mass[atom_type]
            all_weight += atom_type_count * atom_mass / C.AVOGADORO_CONST
        # セル内の密度(g/cm^3)
        density = all_weight / volume
        return density
#--------------------------------------------------------------------------------------
    def count_atom_types(self, res_type='series', condition=None): #k
        """原子のタイプごとに原子の個数をカウントする関数
        Parameters
        ----------
        res_type : str
            res_type='series'のときは結果をpd.Seriesで返す
            res_type='dict'のときは結果をdictで返す
        condition : function
            condition関数
        """
        if condition is None:
            target_atoms = np.array([True] * self.get_total_atoms())
        else:
            target_atoms = condition(self)
        if res_type == 'series':
            return self.atoms.loc[target_atoms,'type'].value_counts().rename(index=self.atom_type_to_symbol)
        elif res_type == 'dict':
            return self.atoms.loc[target_atoms,'type'].value_counts().rename(index=self.atom_type_to_symbol).to_dict()
        else:
            raise ValueError(
                f'res_type: {res_type} is not supported. supported res_type : [series, dict]')
#---------------------------------------------------------------------------------------------
    def shuffle_type(self, type_ratio: list[float]):
        """sfのtypeをランダムにシャッフルする。
            atomsに座標を持たせてから使用。
        Parameters
        ----------
        type_ratio: list[float]
            typeに対する割合が入ったlist

        Example
        -------
        sf.shuffle_type([1,2,3])
            # 原子数:6 -> sf.atoms["type"] = [1,2,2,3,3,3] をシャッフルしたもの
            余りはtype_ratioに応じてランダムに入る
        """
        tot_atoms = self.get_total_atoms()
        tot_ratio = sum(type_ratio)
        type_ratio = [(type_ratio[i]*tot_atoms/tot_ratio) for i in range(len(type_ratio))]
        remain_weight = [type_ratio[i] - int(type_ratio[i]) for i in range(len(type_ratio))]
        type_list = []
        for idx, ratio in enumerate(type_ratio):
            type_list.extend([idx+1 for _ in range(int(ratio))])
        remain_type = random.choices([i+1 for i,_ in enumerate(type_ratio)], k=tot_atoms-len(type_list), weights=remain_weight)
        type_list.extend(remain_type)

        random.shuffle(type_list)    
        self.atoms["type"] = type_list