# limda scripts
ターミナル上で直接実行できるscript
## change_file_type.py
ファイルの形式を変換できる<br>
変換前 : input, car, dumppos, xyz <br>
変換後 : input, dumppos, xyz <br>
```
# dumppos -> input 
change_file_type.py -i dump.pos.5000 -o input.rd
```
-i 変換前のファイル名<br>
-o 変換後のファイル名<br>
-p para, str, limda default対応

## consolidate_dumpposes.py
複数のdumpposを一つのlammps形式のdumpposに変換する.<br><br>
current dirのdumpposを全てまとめてmd.posを作成
```
consolidate_dumpposes.py . 
```
最初に読み込むdirのpathを指定<br>
-o 出力されるlammps dumpposの名前<br>
-p para str, limda_default対応<br>
-s いくつおきにfileを読み込むか<br>

## count_atom_types.py
元素種類ごとに原子数をカウントする.
```
count_atom_types.py -i input.rd
```
-i 読み込むファイル名
-p para str, limda_default対応

## count_bonds.py
結合の種類と数を得る。結合距離はlimda_default["bond_length"]を使用<br>
fileを指定すれば(-i),そのファイルのみ対象、指定しなければcurrent dirのdumppos fileが対象.<br><br>
dump.pos.100のみ対象, dictがprintされる
```
count_bonds.py -i dump.pos.100
```
すべてのdumpposが対象, pandas.DataFrameがprintされる
```
count_bonds.py
```
-i : 読み込むファイル<br>
-p : para<br>
-s : いくつおきにfileを読み込むか<br>

## count_mols.py
分子の数をカウントする.
fileを指定すれば(-i),そのファイルのみ対象、指定しなければcurrent dirのdumppos fileが対象.<br><br>
dump.pos.100のみ対象, dictがprintされる
```
count_mols.py -i dump.pos.100
```
すべてのdumpposが対象, pandas.DataFrameがprintされる
```
count_mols.py
```
-i : 読み込むファイル<br>
-p : para<br>
-s : いくつおきにfileを読み込むか<br>

## density.py
ファイルを読み込み密度を求める.
```
density.py -i dump.pos.100
```
-i : 読み込むファイル<br>
-p : para<br>
