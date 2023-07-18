import random
import pandas as pd
from limda import SimulationFrame

def neighbor_test_case(num: int):
    """
    SimulationFrame の get_neighbor_listのためのtestcaseを作成し、テストします.

    mesh法と愚直なO(N^2)の方法でのneighbor_listを比較します。
    CHOの3種類をランダムな個数(atom_num)、を
    ランダムな座標で(atoms_x,y,z)、
    ランダムな大きさのセル内に入れ(cell[3])、
    ランダムな距離で結合を持つか判定させます(bond_length).
    Condition
    ---------
        cell: list[float] セルサイズ、 10~50 
        mn_cell: float セルの最短の辺の長さ
        atom_num: int 系内の原子数、100~5000 
        atoms_x,y,z list[float] 原子の座標, 0 ~ cell
        bond_length: list[list[float]] 結合の長さ 0 ~ mn_cell-0.01
    """
    t_cnt: int = 0
    for i in range(num):
        sf = SimulationFrame()
        sf.import_para_from_str("C H O")
        cell = [random.uniform(10,50) for _ in range(3)]
        cell_mn = min(cell)
        atom_num = random.randint(100,5000)
        bond_length = [[random.uniform(0, cell_mn/3-0.01) for _ in range(3)] for __ in range(3)]
        bond_length[1][0] = bond_length[0][1]
        bond_length[2][0] = bond_length[0][2]
        bond_length[2][1] = bond_length[1][2]
        atoms_type = [random.choice([1,2,3]) for _ in range(atom_num)]
        atoms_x = [random.uniform(0,cell[0]) for _ in range(atom_num)]       
        atoms_y = [random.uniform(0,cell[1]) for _ in range(atom_num)]       
        atoms_z = [random.uniform(0,cell[2]) for _ in range(atom_num)]
        atoms_list = list(zip(atoms_type, atoms_x, atoms_y, atoms_z))
        sf.atoms = pd.DataFrame(atoms_list, columns=['type', 'x', 'y', 'z'])       
        sf.cell = cell  
        neighbor_list_cy = sf.get_neighbor_list(bond_length = bond_length)
        neighbor_list_py = sf.get_neighbor_list_test(bond_length)
        
        if neighbor_list_cy == neighbor_list_py:
            t_cnt += 1
        else:
            print("-----------------miss---------------------")
            print(f"bond_length\n{bond_length[0]}\n{bond_length[1]}\n{bond_length[2]}")
            print(f"cell | {cell}")
            print(f"atom_num | {atom_num}")
            print("------------------------------------------")

        print(i, flush = True)
    print(f"result : correct ratio is {t_cnt} / {num}")

# 実行         
neighbor_test_case(3300)