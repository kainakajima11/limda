import pandas as pd
import numpy as np

from .SimulationFrame import SimulationFrame

class SimulationFrames(

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

    def __init__(self):
        pass

