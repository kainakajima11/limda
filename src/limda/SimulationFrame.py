import pandas as pd
import numpy as np

from .import_file import ImportFile

class SimulationFrame(
    ImportFile,
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

    def __init__(self):
        self.atoms = None
        self.cell = None
        self.atom_symbol_to_type = None
        self.atom_type_to_symbol = None
        self.atom_type_to_mass = None
        self.step_num = None

    def replicate_atoms(self, replicate_directions=[1, 1, 1]) -> None:
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


