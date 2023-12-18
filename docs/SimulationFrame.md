# SimulationFrame (sf)
sfについての簡単な説明です。<br>
詳しい内容は実際のコードやdocstringを見て下さい。
```python

sf = SimulationFrame() # sfを定義、.limda.yamlにparaがあれば読み込まれる
sf = SimulationFrame("Cr Mn Fe Co Ni") # 定義と同時にparaをimport
sf.import_dumppos("/nfshome17/knakajima/work/dump.pos.100") # 構造を読み込む
```
## 目次
- [変数](#anchor1)<br>
- [メソッド](#anchor2)
  - [SimulationFrame](#anchor3)
  - [ImportFrame](#anchor4)
  - [ExportFrame](#anchor5)
  - [Calculate](#anchor6)
  - [AnalyzeFrame](#anchor7)
<a id="anchor1"></a>
# 変数 
## atoms
原子に関する情報が入ったPandasのDataframe<br><br>
```
>>> sf.atoms
      type       x       y       z     vx     vy     vz    q  mask
0        2  39.410  39.429  -0.324  0.001  0.003  0.001  0.0     0
1        1   0.009   1.876   1.433  0.002 -0.003 -0.001  0.0     0
2        3   1.797  39.322   1.573  0.001  0.003 -0.005  0.0     0
3        5   1.837   1.749  -0.250  0.002  0.001  0.004  0.0     0
4        2   3.650   0.011  -0.294 -0.004 -0.000  0.001  0.0     0
...    ...     ...     ...     ...    ...    ...    ...  ...   ...
9675     4  34.132  37.649  68.430 -0.002  0.002 -0.001  0.0     0
9676     2  35.842  35.869  68.477 -0.004  0.000  0.001  0.0     0
9677     2  35.981  37.684  70.312  0.000  0.002 -0.000  0.0     0
9678     1  37.716  35.817  70.207 -0.000  0.005 -0.001  0.0     0
9679     4  37.678  37.710  68.506  0.001  0.001  0.003  0.0     0

[9680 rows x 9 columns]
```
### Columns
- id : 原子のid
- type : 原子のtype, 1-indexed
- x, y, z : 原子のx, y, z座標
- fx, fy, fz : 原子のx, y, z方向の力
- pred_fx, pred_fy, pred_fz : 原子のNNPによって予測された力
  
## cell
セルのx,y,z方向の大きさが入ったnumpy配列
```
>>> sf.cell
[39.49 39.49 72.518]
```
## atom_symbol_to_type
キーを元素記号、値をtypeとするdict
```
>>> sf.atom_symbol_to_mass
{'Cr': 1, 'Mn': 2, 'Fe': 3, 'Co': 4, 'Ni': 5}
```
## atom_type_to_symbol
キーをtype、値を元素記号とするdict
```
>>> sf.atom_type_to_symbol
{1: 'Cr', 2: 'Mn', 3: 'Fe', 4: 'Co', 5: 'Ni'}
```
## atom_type_to_mass
キーをtype、値を原子の質量とするdict
```
>>> sf.atom_type_to_mass
{1: 51.9961, 2: 54.938045, 3: 55.845, 4: 58.933195, 5: 58.6934}
```
## step_num
mdにおけるstep、int
```
>>> sf.step_num
100
```
## potential_energy
このフレームが持つポテンシャルエネルギー,  float
## pred_potential_energy
NNPによる推論をして得られたポテンシャルエネルギー, 単位はeV, float
## virial_tensor
このフレームの持つvirialテンソル, 単位はeV, shape:(3,3), np.array[float]
## pred_virial_tensor
NNPによる推論をして得られたvirialテンソル, 単位はeV, shape:(3,3), np.array[float]
## limda_default
ImportFrame classが持つインスタンス変数、dict
<a id="anchor2"></a>
# メソッド
<a id="anchor3"></a>
# SimulationFrame 
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/SimulationFrame.py)
## sf[column名]
sf.atoms[column名]と同じ
```
>>> sf["x"] # = sf.atoms["x"]
0       39.410
1        0.009
2        1.797
3        1.837
4        3.650
         ...
9675    34.132
9676    35.842
9677    35.981
9678    37.716
9679    37.678
Name: x, Length: 9680, dtype: float64
```
## len(sf), sf.get_total_atoms()
原子数をintとして得る.
```
>>> len(sf) sf.get_total_atoms()
9680 9680
```
## sf.get_atom_type_set()
フレームにどの原子タイプが含まれているかを得る.
```
>>> sf.get_atom_type_set()
{1, 2, 3, 4, 5} # set[int]
```
## sf.wrap_atoms()
cellからはみ出している原子をcellに入れる。
```python3
sf.wrap_atoms()
```
## sf.replicate_atoms()
sfをx,y,z方向に指定数だけつなげたものに拡張します。
```python3
# x方向に2倍, y方向に3倍, z方向に4倍にする
sf.replicate_atoms([2,3,4]) # [x, y, z]
```
## sf.concat()
2つのsfを結合します。
```python3
sf2 = SimulationFrame("H O") # 別のフレームを用意
sf.concat(sf2) # sfにsf2を結合する
```
## sf.delete_atoms()
条件に当てはまる原子を削除します。
```
>>> condition = lambda i : i["x"] < 30 # x座標が30未満ならTrue
>>> sf.delete_atoms(condition=condition, reindex=False) # 条件式とidxの再配列を行うかを指定
>>> sf["x"]
0       39.410 # x 座標が30未満の原子が消されている
34      30.479 # reindex = False なので idx 1 など欠番になっている
35      30.486
36      32.442
37      32.363
         ...
9675    34.132
9676    35.842
9677    35.981
9678    37.716
9679    37.678
Name: x, Length: 2406, dtype: float64
```
## sf.density()
密度を得る。
```
>>> sf.density()
7.814672036675934 # float
```

## sf.count_atom_types
タイプごとの原子数を得る.
```python3
sf.count_atom_types(res_type="series", condition=condition)
# res_typeは返り値のtype, "series"ならばpandas.Series, "dict"ならdict
# conditionは条件式、指定がなければ全体になる
```
## sf.shuffle_type
タイプをシャッフルする.
```python3
# type1 : type2 : type3 = 1 : 2 : 3 にする.
sf.shuffle_type([1,2,3])
```

## sf.make_magmom_str
vaspのINCARに必要なMAGMOMの値を作成する.
タイプごとにMAGMOMを指定できます.
```python3
# type1,2,3の初期スピンを-1.0, 2.0, 4.0で指定
incar_config["MAGMOM"] = sf.make_magmom_str([-1.0, 2.0, 4.0])
```

## change_lattice_const
格子定数を変える.
```python3
# cellの大きさが[3.0,4.0,5.0]になるように拡大(縮小)
sf.change_lattice_const(new_cell=[3.0,4.0,5.0])

# cellの大きさを2倍に拡大
sf.change_lattice_const(magnification=2.0)
```
new_cell, magnificationのどちらか一方を指定する.

## make_empty_space
真空領域を作成する.
```python3
sf.make_empty_space(empty_length = 10.0   # 真空の大きさは10Å
                    direction = "x"  # "x" or "y" or "z"
                    both_direction = False # Trueのとき両側に真空領域を半分ずつ追加する.
                    )
```
## sf.mirroring_atoms
指定した方向に反転させた構造を結合する.
```python3
sf.mirroring_atoms("x") # "x" or "y" or "z"
# pp -> ppqq みたいな感じ 
```

## sf.shuffle_type_by_part
原子タイプをシャッフルする.<br>
sf.shuffle_typeでは全体で割合が一定になるが,このメソッドでは細かく区切った領域で割合が一定になる.
```python3
sf.shuffle_type_by_part(segment_num = [3,3,3],   # x,y,zを3*3*3個の領域に区切る. 
                        type_ratio = [1,2,3])    # type1:2:3 = 1:2:3にする.
```

## sf.silde_atoms
原子を平行移動させる.
```python3
sf.slide_atoms(slide_length = [1.0,2.0,3.0]  # 原子をx方向に1.0, y方向に2.0, z方向に3.0スライドさせる.
               change_cellsize = True)   # Trueならばスライドした分だけセルサイズも変更する
```
<a id="anchor4"></a>

# ImportFrame
ファイルを読み込んでsfにデータを入れるメソッドが入っています。<br>
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/import_frame.py)

## import_limda_default
~/.limda.yamlを読み込んで,sf.limda_defaultに入れる.<br>
コンストラクタで実行されます.

## import_input()
input fileを読み込みます。<br>
```python3
sf.import_input("/nfshome17/knakajima/work/input.rd")
```

## import_para_from_list()
指定した元素記号からatom_symbol_to_type, atom_type_to_symbol, atom_type_to_massを作成します。<br>
```python3
sf.import_para_from_list(["Cr", "Mn", "Fe", "Co", "Ni"]) # list[str]
```
## import_para_from_str()
機能はimport_para_from_listと同じです。
```python3
sf.import_para_from_list("Cr Mn Fe Co Ni") # 空白区切りのstr
```
## import_car()
car fileを読み込みます。
```python3
sf.import_car("/nfshome17/knakajima/work/Ni.car")
```
## import_dumppos() 
dumppos fileを読み込みます。
```python
sf.import_dumppos("/nfshome17/knakajima/work/dump.pos.100") 
```
## import_mol()
ase ライブラリを用いて、分子を読み込みます。
```python3
sf.import_mol("H2O") # 取り込みたい分子式をstr型で
```
## import_from_poscar()
POSCARからその構造の原子のlist[int]が返される。
```python3
# Crが2つ,Niが4つ入った構造の時
>>> sf.import_from_poscar("/nfshome17/knakajima/work/CrNi/CrNi_0/POSCAR")
[1, 1, 5, 5, 5, 5]    # list[int]
```

## import_xyz()
xyz fileを読み込む
```python3
sf.import_xyz("/nfshome17/knakajima/work/a.xyz")
```

## import_file()
指定した読み込むfileの名前から、fileの種類を分類し読み込む。<br>
読み込めるファイルは
- input file ("input"が含まれる)
- xyz file ("xyz"で終わる)
- car file ("car"で終わる)
- cif file ("cif"で終わる) # TODO
- dump pos file ("dump"か"pos"が含まれる)
```python3
sf.import_file("/nfshome17/knakajima/work/Ni.car") # carで終わっているのでcar fileとしてimportされる.
sf.import_file("/nfshome17/knakajima/work/dump.pos.100") # dumpで始まっているのでdumppos fileとしてimportされる
```

## import_cif
cif fileを読み込む. aseライブラリを用いています.
```python3
sf.import_cif(""/nfshome17/knakajima/work/Cr.cif")
```

<a id="anchor5"></a>
# ExportFrame
主にsfのデータをファイルに出力するメソッドが入っています。<br>
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/export_frame.py)

## export_vasp_poscar()
vaspの計算に必要なPOSCARを作成します。<br>
```python3
sf.export_vasp_poscar(ofn = "POSCAR",    # 出力されるfileのpath
                      comment = "",      # 1行目にかかれるコメント
                      scaling_factor = 1.0)   # 参照vasp wiki
```
## export_vasp_poscar_from_contcar()
vaspの計算に必要なPOSCARを作成します。<br>
sfのデータではなくCONTCARを指定してそれをPOSCARにします。<br>
座標と速度を引き継ぐことができます.
```python3
sf.export_vasp_poscar_from_contcar(ofn = "POSCAR",    # 出力されるfileのpath
                                   contcar_path = "../Ni_0/CONTCAR") # 読み込むCONTCARのpath
```
### export_vasp_incar()
vaspの計算に必要なINCARを作成します。
```python3
sf.export_vasp_incar(ofn="INCAR", # 出力されるfileのpath
                     config=incar_config) # INCARの中身をdictで指定
```
## export_vasp_kpoints()
vaspの計算に必要なKPOINTSを作成します。
```python3
sf.export_vasp_kpoints(ofn="KPOINTS",  # 出力されるfileのpath
                       comment="",     # 1行目にかかれるコメント
                       kx=2,           # x方向のk点
                       ky=2,           # y方向のk点
                       kz=2)           # z方向のk点
```
## export_vasp_iconst()
vaspの計算で使われるICONSTを作成します。
```python3
iconst_config = ['LA 1 2 0', 'LA 1 3 0', 'LA 2 3 0'] # cellを直方体固定
sf.export_vasp_iconst(ofn="ICONST",     # 出力されるfileのpath
                      config=iconst_config) 
```
### export_vasp_potcar()
vaspの計算に必要なPOTCARを作成します。
```python3
sf.export_vasp_poscar(ofn="POSCAR",  # 出力されるfileのpath
                      potcar_root="/nfshome15/POTLIST/potpaw_PBE.54") # potcarのpath
```
## export_dumppos()
sfをdump.pos 形式のファイルとして出力します。
```python3
sf.export_dumppos(ofn="showdump.pos", # showdump.pos fileとして出力
                  time_step=100, # 何ステップ目か
                  out_columns=['type', 'mask', 'x', 'y', 'z']) # 出力されるカラムをlist[str]で指定
```
## export_input()
sfをinput.rd の形式のファイルとしてを出力します。
```python3
sf.export_input(ofn="input.rd", # input.rdで指定
                mask_info = ["#strain - - - - z 1.0"]) # moveなどの情報をlist[str]で指定
```
mask_infoの方式は,laich / lax に従うこと.

## export_xyz()
sfをxyzの形式のファイルとして出力します。
```python3
sf.export_xyz(ofn="a.xyz",
              out_columns=['type', 'x', 'y', 'z'], # 出力されるカラムをlist[str]で指定
              structure="structure") # 構造の名前を指定できる
```
## export_file()
sfを指定したファイル名から形式を判断して、出力します。<br>
出力できるファイルの形式
- input file ("input" で始まる)
- dump.pos file ("dump" or "pos"が含まれる)
- xyz file ("xyz"で終わる)
```python3
sf.export_file("showdump.pos") # dump pos fileとして出力
sf.export_file("input.rd") # input fileとして出力
```
<a id="anchor6"></a>
# Calculate
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/calculate.py)
## vasp()
VASPを実行します.<br>
```python3
sf.vasp(calc_directory="",  # 計算を実行するディレクトリ
        system_name="Cr",     # 計算する系の名前 
        step_num= 2,          # 何回目の実行か # この場合 Cr_2 というディレクトリで計算される
        poscar_comment="",    # poscarの1行目にかかれるコメント
        poscar_scaling_factor=1.0,  # poscarのscaling factor # vasp wiki参照
        incar_config=None,   # incarにかかれる条件をdict[str]で指定
        potcar_root=None,     # potcarのある場所
        kpoints_comment="",   # kpointsの1行目に書かれるコメント
        kpoints_kx = 1, # x方向のk点
        kpoints_ky = 1, # y方向のk点
        kpoints_kz = 1, # z方向のk点
        iconst_config = None, # iconstにかかれる情報, NPTの時は指定する
        vasp_command = "vasp_std", # vaspを実行するためのコマンド
        print_vasp = True,  # 標準出力をpythonからも表示するかどうか
        exist_o = False, # Trueなら同名のディレクトリがあったとき上書きして実行する
        poscar_from_contcar = False, # poscarをcontcarを指定して作成するか
        contcar_path = "", # contcarのpath 
        place = "kbox", # 実行環境 "kbox" or "masamune"
        num_nodes = 1, # 計算ノード数
        )
```
## laich()
Laichを実行します。<br>

## packmol()
Packmolを実行します.<br>
```python
# H2Oを詰めます # セル全体に30個ランダムに詰める.
sf_h2o = SimulationFrame("C H O") 
sf_h2o.import_mol("H2O") 
xyz_condition = [
    [0.7, 0.7, 0.7, sf.cell[0] - 0.7, sf.cell[1] - 0.7, sf.cell[2] - 0.7], # packする範囲を指定 # 端には配置しないようにする 
]
sf.packmol(sf_list=[sf_h2o],
            pack_num_list=[30], # 詰める分子の個数
            tolerance=1.7, # それぞれの水分子は最低1.7Å離れて配置される
            xyz_condition=xyz_condition, # packする範囲条件
            print_packmol=True, # packmolの標準出力を表示するか
            seed=-1) # シード値 # -1 はランダム
```
## lax()
laxを実行します.<br>
```python3
lax_config = {
    "Mode":"MD",
     "NNPModelPath":"/nfshome15/rkudo/work/lax_bench/train/model1/frozen_models/allegro_frozen_544000.pth",
     "TotalStep":"100000",
     "TimeStep":"0.50",
     "FileStep":"100",
     "SaveRestartStep":"0",
     "MPIGridX":"1",
     "MPIGridY":"1",
     "MPIGridZ":"1",
     "CUTOFF":"4.0",
     "MARGIN":"1.0",
     "GhostFactor":"20.0",
     "InitTemp":"300.0",
     "Thermo":"NoseHoover",
     "AimTemp":"300.0",
     "FinalTemp":"300.0",
     "ThermoFreq":"10.0",
     "Baro":"NoseHoover",
     "PressFlag":"1 1 1",
     "AimPress":"1.0 1.0 1.0",
     "BaroFreq":"3000 3000 3000",
     "FlagArrayUsage":"1",
     "MaxPairNumberOfOneMeshDevided":"3000",
}
sf.lax(calc_dir = "lax_calc/lax_calc_0", # 計算が行われるdir
       lax_cmd = "~/lax/src/build/lax", # 実行可能なlaxのpath
       lax_config = lax_config, # config.rdに必要な内容をdictで渡す.
       print_lax = True,  # out fileを出力するか
       exist_ok = False, # calc_dirが存在するときに実行するか
       mask_info = [] # input.rdに書かれる#pressz、#moveなどの情報
       # ex. mask_info = [#pressz 1 1 0", "#move 2 x 100 - - - -"] # mask 1にpress, mask 2にmove
       omp_num_threads = 3 # OMP_NUM_THREADSの値
       )

```
## allegro()
sfに入っている原子の座標などから、かかる力とポテンシャルエネルギーを予測します.<br>
力はatomsの["pred_x", "pred_y", "pred_z"]に入ります。<br>
ポテンシャルエネルギーはpred_potential_energyに入ります。<br>
Virialテンソルの値はpred_virial_tensorに入ります。shapeは(3,3) <br>
```python3
sf.allegro(
    cut_off = 4.0, #cutoff
    device = "cuda", #device, 'cpu' or 'cuda'
    allegro_model = torch.jit.load("allegro_frozen_800000.pth")
    # frozenされたAllegroを読み込んだモデル # pathではない
    flag_calc_virial = False, # virialを推論するか
    )
```
<a id="anchor7"></a>
# AnalyzeFrame
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/analyse_frame.py)

## get_neighbor_list
sfに対する、隣接リストを作成し、返します。<br>
neighor_listの形式はlist[list[int]です.<br>
配列のi番目の要素はi番目の原子と隣接する原子のidxの配列です.<br>
mode = "bond_length"ならば元素種類ごとに長さを指定する. <br>
mode = "cut_off"ならばfloat型の一つの値で長さを定義.
```python3
# bond_length # 2元素系
neighbor_list = sf.get_neighbor_list(mode="bond_length", bond_length=[[1.2, 2.0],[2.0, 2.3]])
# cut_off
neighbor_list = sf.get_neighbor_list(mode="cut_off", cut_off=3.4)
```
## get_edge_idx
隣接リストをallegroのデータセットの形式にしたもの(edge_idx)を返します。<br>
edge_idxはlist[list[int, int]]で、配列の要素は大きさ2の配列であり、<br>
[a,b]のときa番目の原子とb番目の原子は隣接してることを表します。
```python3
edge_index = sf.get_egde_index(cut_off=3.4)
```
## get_mols_list
分子ごとに原子のidをlistにまとめる. <br>
例えば、水分子が3個とアンモニアが1個の系では
[[0, 1, 2],  # 水分子 <br>
  [3, 4, 5],  # 水分子 <br>
  [6, 7, 8],  # 水分子 <br>
  [9, 10, 11, 12]] # アンモニア <br>
のようなlist[list[int]]が得られる.
```python3
mols_list = sf.get_mols_list(mode="bond_length", bond_length=[[1.2, 2.0],[2.0, 2.3]])
```

## get_mols_dict
分子種類ごと, 分子ごとに原子のidをまとめる. <br>
例えば、水分子が3個とアンモニアが1個あるときは <br>
    {"H2O1":[[0, 1, 2], [3, 4, 5], [6, 7, 8]], <br>
     "H3N1":[[9, 10, 11, 12]]} というdict[str, list[list[int]]]が得られる.  <br>
```python3
mols_dict = sf.get_mols_dict(mode="bond_length", bond_length=[[1.2, 2.0],[2.0, 2.3]])
```
## count_mols
sf内に何の分子が何個存在するかを返します。<br>
dict[str, int]の形式で返され
H2O1 : 3 であれば,水分子が3個あることを表します。<br>
```python3
mols_count = sf.count_mols(mode="bond_length", bond_length=[[1.2, 2.0],[2.0, 2.3]])
```
## count_bonds
sf内に何の結合が何個あるかを返します。<br>
dict[str, int]の形式で返され、
H-O : 2 であれば水素酸素結合が2つ存在することを表します。
```python3
count_bonds_dict = sf.count_bonds(mode="bond_length", bond_length=[[1.2, 2.0],[2.0, 2.3]])
```

## get_sum_of_momentums
各方向の運動量の合計を計算する.
[x,y,z]の運動量がの合計が入ったndarrayが得られる.
```python3
momentum_sum = sf.get_sum_of_momentums()
```
