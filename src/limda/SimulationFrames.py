import pandas as pd
import numpy as np
import pathlib
import random
import os
from tqdm import tqdm, trange
from typing import Union
import torch
from .import_frames import ImportFrames
from .export_frames import ExportFrames
from .SimulationFrame import SimulationFrame
from . import Default as D

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
    def __init__(self, para: str=""):
        self.sf:list[SimulationFrame] = []
        self.atom_symbol_to_type: dict[str, int] = None
        self.atom_type_to_symbol : dict[int, str] = None
        self.atom_type_to_mass : dict[int, float] = None

        D.set_default()

        if para:
            self.import_para_from_str(para)
        elif D.PARA is not None:
            self.import_para_from_list(D.PARA)
#---------------------
    def __len__(self):
        """
        len(sfs)でlen(sfs.sf)を得ることができる
        """
        return len(self.sf)
#----------------------------- 
    def __getitem__(self, key):
        """sfs = SimulationFrames()
        sfs[step_idx]でsfs.sf[step_idx]を得ることができる
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
        
        for step_idx, frame in enumerate(self.sf):
            frame.step_num = step_idx
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

#---------------------------------------------------------------------------------------------------------------
    def allegro(self,
                cut_off: float,
                device: Union[str, torch.DeviceObjType],
                allegro_model: torch.jit._script.RecursiveScriptModule,
                flag_calc_virial:bool = False,
                ) -> None:
        """sfsのもつすべてのSimulationFrameに対して、以下を行う
        Allegroを使って、sfに入っている原子の座標に対して推論を行い、
        ポテンシャルエネルギーと原子に働く力を計算する
        Allegroによって予測されたポテンシャルエネルギーはsf.pred_potential_energyに入る
        Allegroによって予測されたそれぞれの原子にかかる力はsf.atoms.loc[:, ["pred_fx", "pred_fy", "pred_fz"]]に入る
        Parameters
        ----------
            cut_off: float
                原子の相互作用のカットオフ半径
            device: Union[str, torch.DeviceObjType]
                どのデバイスでAllegroの計算を行うか, "cpu" or "cuda"
            allegro_model: torch.jit._script.RecursiveScriptModule
                frozenされたAllegroを読み込んだモデル
                pathではないことに注意
        """
        for frame_idx in range(len(self.sf)):
            self.sf[frame_idx].allegro(cut_off=cut_off,
                                       device=device,
                                       allegro_model=allegro_model,
                                       flag_calc_virial=flag_calc_virial)

#---------------------------------------------------------------------------------------------- 
    def concat_force_and_pred_force(self, 
                                    reduce_direction: bool = False,
                                    ) -> pd.DataFrame:
        """それぞれのSimulationFrameが持つ力とAllegroによって予測された力という一つのDataFrameを作る
        Parameters
        ----------
            reduce_direction: bool = False
                x, y, z方向をまとめるかどうか
                reduce_direction == Trueのときは、返り値のDataFrameのcolumnは["f", "pred_f"]

                reduce_direction == Falseのときは、返り値のDataFrameのcolumnは
                ["fx", "fy", "fz", "pred_fx", "pred_fy", "pred_fz"]
        Returns
        -------
            force_and_pred_force_reduced または force_and_pred_forceが返る
            force_and_pred_force_reduced : pd.DataFrame
                columnは["f", "pred_f"]、行は同じ原子の同じ方向に働く力に対応する
            force_and_pred_force : pd.DataFrame
                columnは["fx", "fy", "fz", "pred_fx", "pred_fy", "pred_fz"]、行は同じ原子の同じ方向に働く力に対応する

        Usage
        -----
            横軸に第一原理計算で計算された力, 縦軸にAllegroで予測された力というプロットと作りたいときに用いる
            ```python
            sfs = SimulationFrames("C H O N Si")
            sfs.import_allegro_frames("/path/to/allegro_frame.pickle")
            device = torch.device("cuda")
            model = torch.jit.load("/path/to/allegro_frozen_1000.pth")
            model.to(device)
            sfs.allegro(cut_off=3.4, allegro_model=model, device=device)
            force_and_pred_force = sfs.concat_force_and_pred_force(reduce_direction=True)
            plt.scatter(force_and_pred_force.loc[:, "f"], force_and_pred_force.loc[:, "pred_f"], alpha=0.005)
            plt.show()
            ```
        """
        force_and_pred_force_list = []
        for frame_idx in range(len(self.sf)):
            force_and_pred_force_list.append(
                self.sf[frame_idx].atoms.loc[:,
                    ["fx", "fy", "fz", "pred_fx", "pred_fy", "pred_fz"]
                ]
            )
        force_and_pred_force = pd.concat(force_and_pred_force_list)

        if reduce_direction:
            force_and_pred_force_reduced = pd.DataFrame(
                np.concatenate(
                    [force_and_pred_force.loc[:, ["fx", "pred_fx"]].values,
                    force_and_pred_force.loc[:, ["fy", "pred_fy"]].values,
                    force_and_pred_force.loc[:, ["fz", "pred_fz"]].values,],
                    axis=0),
                columns=["f", "pred_f"]
            )
            return force_and_pred_force_reduced
        else:
            return force_and_pred_force
#----------------------------------------------------------------------------------------------
    def concat_pot_and_pred_pot(self) -> pd.DataFrame:
        """それぞれのSimulationFrameが持つポテンシャルエネルギーと
        Allegroによって予測されたポテンシャルエネルギーという一つのDataFrameを作る
        Returns
        -------
            pot_and_pred_pot : pd.DataFrame
                columnは["potential_energy", "pred_potential_energy"]、行は同じフレームに対応する
        Usage
        -----
            横軸に第一原理計算で計算されたポテンシャルエネルギー, 
            縦軸にAllegroで予測されたポテンシャルエネルギーというプロットと作りたいときに用いる
        """
        potential_energy_list = []
        pred_potential_energy_list = []
        for frame_idx in range(len(self.sf)):
            potential_energy_list.append(self.sf[frame_idx].potential_energy)
            pred_potential_energy_list.append(self.sf[frame_idx].pred_potential_energy)

        pot_and_pred_pot = pd.DataFrame({
                "potential_energy":potential_energy_list,
                "pred_potential_energy":pred_potential_energy_list
            })

        return pot_and_pred_pot