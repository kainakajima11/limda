# limda default
limdaでよく使う変数などをデフォルトとして使える機能<br>
~/.limda.yamlというyaml fileを各自で作成してください。
```python3
sf = SimulationFrame()
sfs = SimulationFrames()
```
sf、もしくはsfsを定義したときに自動で.limda.yamlが読み込まれます。<br>
yamlファイルの中身は インスタンス変数limda_defaultにdictとして入ります。

## メソッド内で使用されるkey
### para
paraを表すlist[str], paraを読み込む際に引数を空にすればデフォルトが使用されます。
```yaml
para : ["Cr", "Mn", "Fe", "Co", "Ni"]
```
defaultがある場合、以下の二つのコードは同じ働きになります。
```python3
sf = SimulationFrame()
```
```python3
sf = SimulationFrame()
sf.import_para_from_list(["Cr", "Mn", "Fe", "Co", "Ni"]
```
### cut_off
cut_offの値, float<br>
neighbor_listをcut_offで作成するときに使われます。
```yaml
cut_off : 4.0
```
defaultがある場合、以下のコードは同じ働きになります。
```python3
edge_index = sf.get_edge_index()
```
```python3
edge_index = sf.get_edge_index(cut_off=4.0)
```
### bond_length
元素種類ごとの結合の長さを表すlist[list[float]]<br>
neighbor_listをbond_lengthで作成するときに使われます。
```yaml
bond_length :
 - [1.1, 1.2, 1.3, 1.4, 1.5]
 - [1.2, 2.2, 2.3, 2.4, 2.5]
 - [1.3, 2.3, 3.3, 3.4, 3.5]
 - [1.4, 2.4, 3.4, 4.4, 4.5]
 - [1.5, 2.5, 3.5, 4.5, 5.5]
```
defaultがあル場合、以下のコードで結合の種類と数が得られます。
```python3
bonds_dict = sf.count_bonds(mode="bond_length")
```
### NELM
vaspで1stepあたりの最大iteration回数
```yaml
NELM : 400
```
defaultがある場合、datasetを作成するときに収束しなかった構造を勝手に除きます。

## その他の活用
メソッド以外でよく使う変数を入れておく活用ができます。<br>
例: kboxとmasamuneで異なる変数をdefault化
```yaml
# kbox
car_dir_path : "/nfshome17/knakajima/CAR"
# masamune
car_dir_path : "/home/kainkjm/CAR"
```
```python3
sf.import_car(pathlib.Path(sf.limda_default["car_dir_path"]) / "Ni.car")
```
どちらの環境でも同様に機能するようになります。

##### TODO vasp_incar_config


