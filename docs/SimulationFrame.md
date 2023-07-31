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
Columns
- id : 原子のid
- type : 原子のtype, 1-indexed
- x, y, z : 原子のx, y, z座標
- fx, fy, fz : 原子のx, y, z方向の速度 
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
<a id="anchor2"></a>
# メソッド
<a id="anchor3"></a>
## SimulationFrame 
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/SimulationFrame.py)

### sf[column名]
sf.atoms[column名]と同じ

### len(sf), sf.get_total_atoms()
sfにある原子の総数を返り値として得られる。

### sf.get_atom_type_set()
sfの原子のtypeのsetを返り値として得られる。

### sf.wrap_atoms()
cellからはみ出している原子をcellに入れる。

### sf.replicate_atoms()
引数：x,y,z方向に何倍するかを表すlist<br><br>
sfをx,y,z方向に指定数だけつなげたものに拡張します。

### sf.concat()
引数：sfに結合するSimulationFrame<br><br>
2つのsfを結合します。

### sf.delete_atoms()
引数：
- 削除する原子の条件,関数
- 削除後idを付け直すか,bool
<br>
条件に当てはまる原子を削除します。

### sf.density()
引数：x,y,zの範囲の下限と上限,それぞれfloat
sfの密度を返します。
<a id="anchor4"></a>
## ImportFrame
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/import_frame.py)
<a id="anchor5"></a>
## ExportFrame
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/export_frame.py)
<a id="anchor6"></a>
## Calculate
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/calculate.py)
<a id="anchor7"></a>
## AnalyzeFrame
[実際のコード](https://github.com/kainakajima11/limda/blob/main/src/limda/analyse_frame.py)
