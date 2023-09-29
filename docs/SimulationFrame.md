# SimulationFrame (sf)
sfについての簡単な説明です。<br>
詳しい内容は実際のコードやdocstringを見て下さい。
```python
sf = SimulationFrame() # ただsfをSimulationFrameとして定義
sf = SimulationFrame("Cr Mn Fe Co Ni") # 定義と同時にsf.import_para_from_str()を実行
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
### atoms
原子に関する情報が入ったPandasのDataframe<br><br>
```python3
print(sf.atoms)
```
```
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
#### Columns
- id : 原子のid
- type : 原子のtype, 1-indexed
- x, y, z : 原子のx, y, z座標
- fx, fy, fz : 原子のx, y, z方向の力
- pred_fx, pred_fy, pred_fz : 原子のNNPによって予測された力
### cell
セルのx,y,z方向の大きさが入ったnumpy配列
```python3
print(sf.cell)
```
```
[39.49 39.49 72.518]
```
### atom_symbol_to_type
キーを元素記号、値をtypeとするdict
```python3
print(sf.atom_symbol_to_mass)
```
```
{'Cr': 1, 'Mn': 2, 'Fe': 3, 'Co': 4, 'Ni': 5}
```
### atom_type_to_symbol
キーをtype、値を元素記号とするdict
```python3
print(sf.atom_type_to_symbol)
```
```
{1: 'Cr', 2: 'Mn', 3: 'Fe', 4: 'Co', 5: 'Ni'}
```
### atom_type_to_mass
キーをtype、値を原子の質量とするdict
```python3
print(sf.atom_type_to_mass)
```
```
{1: 51.9961, 2: 54.938045, 3: 55.845, 4: 58.933195, 5: 58.6934}
```
### step_num
laichなどで扱うstep数、int
### potential_energy
sfのポテンシャルエネルギー、float
### pred_potential_energy
sfをNNPで推論して得られたエネルギー、float
<a id="anchor2"></a>
# メソッド
<a id="anchor3"></a>
## SimulationFrame 
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/SimulationFrame.py)
引数に元素記号のstrを入れれば、import_para_from_str()が呼び出されます.

### sf[column名]
sf.atoms[column名]と同じ
```python3
print(sf["x"])
```
```
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
### len(sf), sf.get_total_atoms()
sfにある原子の総数を返り値として得られる。
```python3
print(len(sf), sf.get_total_atoms())
```
```
9680 9680
```
### sf.get_atom_type_set()
sfの原子のtypeのsetを返り値として得られる。
```python3
print(sf.get_atom_type_set())
```
```
{1, 2, 3, 4, 5}
```
### sf.wrap_atoms()
cellからはみ出している原子をcellに入れる。
```python3
sf.cell[0] = 30 # cellのx方向を縮めて、セルから原子がはみ出ている状態にします.
sf.wrap_atoms()
print(sf["x"])
```
```
0        9.41 # 元々39.41だった原子がセルの中に入っている.
1       0.009
2       1.797
3       1.837
4        3.65
        ...
9675    4.132
9676    5.842
9677    5.981
9678    7.716
9679    7.678
```
### sf.replicate_atoms()
sfをx,y,z方向に指定数だけつなげたものに拡張します。
```python3
sf.replicate_atoms([2,2,2]) # [x方向への拡張倍率, y方向への拡張倍率, z方向への拡張倍率]
print(len(sf))
```
```
77440 # 9680 * 2 * 2 * 2
```
### sf.concat()
2つのsfを結合します。
```python3
sf2 = SimulationFrame("H O")
sf.concat(sf2)
```
### sf.delete_atoms()
条件に当てはまる原子を削除します。
```python3
dlt = lambda i : i["x"] < 30 # x座標が30未満ならTrue
sf.delete_atoms(condition=dlt, reindex=False) # 条件式とidxの再配列を行うかを指定
print(sf["x"])
```
```
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
### sf.density()
sfの密度を返します。
```python3
print(sf.density())
```
```
7.814672036675934
```
<a id="anchor4"></a>
## ImportFrame
主にファイルを読み込んでsfにデータを入れるメソッドが入っています。<br>
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/import_frame.py)

### import_input()
input fileを読み込みます。<br>
```python3
sf.import_input("/nfshome17/knakajima/work/input.rd")
```

### import_para_from_list()
指定した元素記号からatom_symbol_to_type, atom_type_to_symbol, atom_type_to_massを作成します。<br>
['C', 'H', 'O']のように引数はlist
```python3
sf.import_para_from_list(["Cr", "Mn", "Fe", "Co", "Ni"])
```
### import_para_from_str()
機能はimport_para_from_listと同じです。
"C H O"のように引数は空白区切りのstr
```python3
sf.import_para_from_list("Cr Mn Fe Co Ni")
```
### import_car()
car fileを読み込みます。
```python3
sf.import_car("/nfshome17/knakajima/work/Ni.car")
```
### import_dumppos() 
dumppos fileを読み込みます。
```python
sf.import_dumppos("/nfshome17/knakajima/work/dump.pos.100") 
```
### import_mol()
ase ライブラリを用いて、分子を読み込みます。
```python3
sf.import_mol("H2O") # 取り込みたい分子式をstr型で
```
### import_from_poscar()
POSCARからその構造の原子のlist[int]が返される。
```python3
# Crが2つ,Niが4つ入った構造の時
print(sf.import_from_poscar("/nfshome17/knakajima/work/CrNi/CrNi_0/POSCAR"))
```
```
[1, 1, 5, 5, 5, 5]
```

### import_xyz()
xyz fileを読み込む
```python3
sf.import_xyz("/nfshome17/knakajima/work/a.xyz")
```

### import_file()
指定した読み込むfileの名前から、fileの種類を分類し読み込む。<br>
読み込めるファイルは
- input file ("input"から始まる)
- dump pos file ("dump"で始まるか"pos"で終わる)
- xyz file ("xyz"で終わる)
- car file ("car"で終わる)
```python3
sf.import_file("/nfshome17/knakajima/work/Ni.car") # carで終わっているのでcar fileとしてimportされる.
sf.import_file("/nfshome17/knakajima/work/dump.pos.100") # dumpで始まっているのでdumppos fileとしてimportされる
```
<a id="anchor5"></a>
## ExportFrame
主にsfのデータをファイルに出力するメソッドが入っています。<br>
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/export_frame.py)

### export_vasp_poscar()
vaspの計算に必要なPOSCARを作成します。<br>
sfのデータが使われます。

### export_vasp_poscar_from_contcar()
vaspの計算に必要なPOSCARを作成します。<br>
sfのデータではなく、CONTCARを指定して、それをPOSCARにします。<br>

### export_vasp_incar()
vaspの計算に必要なINCARを作成します。<br>

### export_vasp_kpoints()
vaspの計算に必要なKPOINTSを作成します。

### export_vasp_iconst()
vaspの計算で使われるICONSTを作成します。

### export_vasp_potcar()
vaspの計算に必要なPOTCARを作成します。

### export_dumppos()
sfをdump.pos 形式のファイルとして出力します。

### export_input()
sfをinput.rd の形式のファイルとしてを出力します。

### export_xyz()
sfをxyzの形式のファイルとして出力します。

### export_file()
sfを指定したファイル名から形式を判断して、出力します。<br>
出力できるファイルの形式
- input file ("input" で始まる)
- dump.pos file ("dump"で始まる or "pos"で終わる)
- xyz file ("xyz"で終わる)

<a id="anchor6"></a>
## Calculate
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/calculate.py)
### vasp()
VASPの計算を実行します.<br>

### laich()
Laichの計算を実行します。<br>

### packmol()
Packmolの計算を実行します.<br>
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
### lax()
laxの計算を実行します.<br>

### allegro()
sfに入っている原子の座標などから、かかる力とポテンシャルエネルギーを予測します.<br>
力はatomsの["pred_x", "pred_y", "pred_z"]に入ります。<br>
ポテンシャルエネルギーはpred_potential_energyに入ります。<br>

<a id="anchor7"></a>
## AnalyzeFrame
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/analyse_frame.py)

### get_neighbor_list()
sfに対する、隣接リストを作成し、返します。<br>
neighor_listの形式はlist[list[int]です.<br>
配列のi番目の要素はi番目の原子と隣接する原子のidxの配列です.<br>
原子のタイプごとに長さを指定して結合リストを作ることもできます.<br>

### get_edge_idx()
隣接リストをallegroのデータセットの形式にしたもの(edge_idx)を返します。<br>
edge_idxはlist[list[int, int]]で、配列の要素は大きさ2の配列であり、<br>
[a,b]のときa番目の原子とb番目の原子は隣接してることを表します。

### count_molecules()
sf内に何の分子が何個存在するかを返します。<br>
dict[str, int]の形式で返され
H2O1 : 3 であれば,水分子が3個あることを表します。<br>

### count_bonds()
sf内に何の結合が何個あるかを返します。<br>
dict[str, int]の形式で返され、
H-O : 2 であれば水素酸素結合が2つ存在することを表します。

### count_coord_numbers()
sf内の原子の配位数の数を返します。
list[dict[int,int]]の形式で返され、i番目の配列が、{1 : 2, 2 : 3}であれば、<br>
i番目のタイプの原子は1配位が2つ、2配位が3つあることを表します.<br>
