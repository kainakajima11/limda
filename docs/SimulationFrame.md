# SimulationFrame (sf)
sfについての簡単な説明です。<br>
詳しい内容は実際のコードやdocstringを見て下さい。

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
#### Columns
- id : 原子のid
- type : 原子のtype, 1-indexed
- x, y, z : 原子のx, y, z座標
- fx, fy, fz : 原子のx, y, z方向の力
- pred_fx, pred_fy, pred_fz : 原子のNNPによって予測された力
### cell
セルのx,y,z方向の大きさが入ったnumpy配列
### atom_symbol_to_type
キーを元素記号、値をtypeとするdict
### atom_type_to_symbol
キーをtype、値を元素記号とするdict
### atom_type_to_mass
キーをtype、値を原子の質量とするdict
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

### len(sf), sf.get_total_atoms()
sfにある原子の総数を返り値として得られる。

### sf.get_atom_type_set()
sfの原子のtypeのsetを返り値として得られる。

### sf.wrap_atoms()
cellからはみ出している原子をcellに入れる。

### sf.replicate_atoms()
sfをx,y,z方向に指定数だけつなげたものに拡張します。

### sf.concat()
2つのsfを結合します。

### sf.delete_atoms()
条件に当てはまる原子を削除します。

### sf.density()
sfの密度を返します。

<a id="anchor4"></a>
## ImportFrame
主にファイルを読み込んでsfにデータを入れるメソッドが入っています。<br>
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/import_frame.py)

### import_input()
input fileを読み込みます。<br>

### import_para_from_list()
指定した元素記号からatom_symbol_to_type, atom_type_to_symbol, atom_type_to_massを作成します。<br>
['C', 'H', 'O']のように引数はlist
### import_para_from_str()
機能はimport_para_from_listと同じです。
"C H O"のように引数は空白区切りのstr
### import_car()
car fileを読み込みます。

### import_dumppos()
dumppos fileを読み込みます。

### import_mol()
ase ライブラリを用いて、分子を読み込みます。

### import_from_poscar()
POSCARから原子の種類と個数がlist[int]で返される。

### import_xyz()
xyz fileを読み込む

### import_file()
指定した読み込むfileの名前から、fileの種類を分類し読み込む。<br>
読み込めるファイルは
- input file ("input"から始まる)
- dump pos file ("dump"で始まるか"pos"で終わる)
- xyz file ("xyz"で終わる)
- car file ("car"で終わる)
  
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

### lax()
laxの計算を実行します.<br>

### allegro()
sfに入っている原子の座標などから、かかる力とポテンシャルエネルギーを予測します.<br>
力はatomsの["pred_x", "pred_y", "pred_z"]に入ります。<br>
ポテンシャルエネルギーはpred_potential_energyに入ります。<br>

<a id="anchor7"></a>
## AnalyzeFrame
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/analyse_frame.py)
