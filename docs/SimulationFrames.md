# SimulationFrames (sfs)
sfsについての簡単な説明です。<br>
詳しい内容は実際のコードやdocstringを見て下さい。
```python3
from limda import SimulationFrames # import
sfs = SimulationFrames() # define sfs
sfs.import_dumpposes("/nfshome17/knakajima/sfstest") # 複数のframeをinput
# /nfshome17/knakajima/sfstest には dump.pos.0 dump.pos.100 dump.pos.200 がある.
```
## 目次
- [変数](#anchor1)<br>
- [メソッド](#anchor2)
  - [SimulationFrame](#anchor3)
  - [ImportFrames](#anchor4)
  - [ExportFrames](#anchor5)
  - [AnalyzeFrames](#anchor6)
<a id="anchor1"></a>
# 変数 
## sf
sfsが何個のフレームを持つか, SimulationFrameのlist
## atom_symbol_to_type
キーを元素記号、値をtypeとするdict
```
>>> sfs.atom_symbol_to_type
{'Cr': 1, 'Mn': 2, 'Fe': 3, 'Co': 4, 'Ni': 5, 'H': 6, 'O': 7}
```
## atom_type_to_symbol
キーをtype、値を原子の質量とするdict
```
>>> sfs.atom_type_to_symbol
{1: 'Cr', 2: 'Mn', 3: 'Fe', 4: 'Co', 5: 'Ni', 6: 'H', 7: 'O'}
```
## atom_type_to_mass
キーをtype、値を原子の質量とするdict
```
>>> sfs.atom_type_to_mass
{1: 51.9961, 2: 54.938045, 3: 55.845, 4: 58.933195, 5: 58.6934, 6: 1.00794, 7: 15.9994}
```
## limda_default
ImportFramesのインスタンス変数,~/.limda.yamlの中身をdictとして持つ

<a id="anchor2"></a>
# メソッド
<a id="anchor3"></a>
## SimulationFrames 
[実際のコード](https://github.com/kainakajima11/limda/src/limda/SimulationFrames.py)
### len(sfs)
len(sfs.sf)と同じ、SimulationFramesが何個のフレームから成るか<br>

## sfs[step数]
sfs.sf[step数]と同じ<br>

## shuffle_sfs
sfs.sfをシャッフルする。<br>
```python3
for frame_i in sfs.sf:
    print(frame_i.step_num) # 0 100 200
sfs.shuffle_sfs()
for frame_i in sfs.sf:
    print(frame_i.step_num) # 100 200 0
```
## concat_sfs
sfsのlistを引数とし、それらを結合する.
```python3
sfs = SimulationFrames()
sfs_1 = SimulationFrames()
sfs_1.import_dumpposes("/nfshome17/knakajima/work/md/npt_test/sfstest")
sfs_2 = SimulationFrames()
sfs_2.import_dumpposes("/nfshome17/knakajima/work/md/npt_test/sfstest")
print(len(sfs_1), len(sfs_2)) # 3 3
sfs.concat_sfs([sfs_1, sfs_2])
print(len(sfs)) # 6
```
## split_sfs_specified_list_size
sfsを複数のsfsに分割する。何分割するかを指定.<br>
返り値として分割されたsfsのリストが得られる。<br>
```python3
# len(sfs) == 6 のとき
sfs_list = sfs.split_sfs_specified_list_size(3) # 3つのsfsに分割
print(len(sfs_list[0]), len(sfs_list[1]), len(sfs_list[2])) # 2 2 2 
```
## split_sfs
sfsを複数のsfsに分割する。一つのsfs.sfの大きさを指定。<br>
返り値として分割されたsfsのリストが得られる。<br>
```python3
# len(sfs) == 6
sfs_list = sfs.split_sfs(3) # 1つのsfsの大きさが3になるように分割
print(len(sfs_list[0]), len(sfs_list[1])) # 3 3 
```
## allegro
allegroを用いて、sfs.sfそれぞれのポテンシャルエネルギー･力･virialテンソルを推論する.<br>
sfs.sf[].pred_potential_energyにポテンシャルエネルギーが、sfs.sf[].loc[:, ["pred_fx", "pred_fy", "pred_fz"]]に力が、sfs.sf[].pred_virial_tensorにvirialテンソルが入る。
```python3
sfs.allegro(
    cut_off = 4.0, #cutoff
    device = "cuda", #device, 'cpu' or 'cuda'
    allegro_model = torch.jit.load("allegro_frozen_800000.pth")
    # frozenされたAllegroを読み込んだモデル # pathではない
    flag_calc_virial = False, # virialを推論するか
    )
```
## concat_force_and_pred_force
力(sfs.sf[].loc[:, ["fx", "fy", "fz"]])と推論された力(sfs.sf[].loc[:, ["pred_fx", "pred_fy", "pred_fz"]])をまとめたDataFrameを作成する.
```python3
sfs.concat_force_and_pred_force(reduce_direction=True)
# reduce_direction がFalseならば, colomnは["fx", "fy", "fz", "pred_fx", "pred_fy", "pred_fz"]
# Trueならばそれらを結合してcolumnは["f", "pred_f"]になる.
```
## concat_pot_and_pred_pot
ポテンシャルエネルギー(sfs.sf[].potential_energy)と推論されたポテンシャルエネルギー(sfs.sf[].pred_potential_energy)をまとめたDataFrameを作成する.<br>
columnは["potential_energy", "pred_potential_energy"]
```python3
sfs.concat_pot_and_pred_pot() 
```
## concat_virial_and_pred_virial
virialテンソル(sfs.sf[].virial_tensor)と推論されたvirialテンソル(sfs.sf[].virial_tensor)をまとめたDataFrameを作成する.<br>
単位がeVからeV/Å^3になることに注意
columnは["stress", "pred_stress"]
```python3
sfs.concar_virial_and_pred_virial(only_diag=False) # Trueならば非対角成分を捨てる
```
## get_step_nums
ステップ数のリストを作成する.
```python3
step_list = sfs.get_step_nums()
```

<a id="anchor4"></a>

# ImportFrames
[実際のコード](https://github.com/kainakajima11/limda/src/limda/import_frames.py)

## import_limda_default
~/.limda.yamlを読み込み、中身をsfs.limda_defaultに格納する.
```python3
sfs.import_limda_default()
```

## import_vasp
vaspで計算した第一原理MDファイルから、原子の座標, cellの大きさ, 原子にかかる力, ポテンシャルエネルギー, virialテンソルを読み込む.
```python3
sfs.import_vasp(calc_directory="Cr_0", # vaspの計算ディレクトリのパス
                NELM=400) # vaspの最大Iteration回数
# NELMは,指定しなければlimda_defaultの値が使用されます.
# NELM回まで達した構造は収束していないと判断し,importしません.
```
## import_dumpposes
フォルダのpathを指定し、そこにあるdumppos fileを読み込む.
```python3
sfs.import_dumpposes(dir_name="/nfshome17/knakajima/work/MD_Cr", # dumpposの入っているフォルダ
                     skip_num=10) # いくつおきdumpposを読み込むのか
```

## import_para_from_list
listからparaを読み込みatom_symbol_to_type, atom_type_to_symbol, atom_type_to_massを作成する.
```python3
sfs.import_para_from_list(["Cr", "Mn", "Fe", "Co", "Ni"]) # list[str]
```
## import_para_from_str
strからparaを読み込みatom_symbol_to_type, atom_type_to_symbol, atom_type_to_massを作成する.
機能はimport_para_from_listと同じです。
```python3
sfs.import_para_from_str("Cr Mn Fe Co Ni") # 空白区切りのstr
```
## import_allegro_frames
pickle fileからcell, position, force, atom_types, potential_energyを得る.<br>
edge_indexとcut_offの情報は抜け落ちる.
```python3
sfs.import_allegro_frames("Cr_0.pickle")
```

<a id="anchor5"></a>

# ExportFrames
[実際のコード](https://github.com/kainakajima11/limda/src/limda/export_frames.py)

## export_dumpposes
sfs.sfの持つ構造それぞれをdumppos fileに出力する.
```python3
sfs.export_dumpposes(output_folder="/nfshome17/knakajima/work/MD_Cr" # 出力するフォルダ
                     out_columns=['type', 'mask', 'x', 'y', 'z'])    # 出力されるcolumn
```

## export_allegro_frames
sfs.sfの持つ構造をpickle fileに出力する.出力されるpickle fileはcell, position, cutoff, edge_index, type, force, potential_energy, virial_tensorの情報を持つ.
```python3
sfs.export_allegro_frames(output_dir="/nfshome17/knakajima/work/dataset/", # 出力するパス
                          output_file_name="Cr", # 出力されるfile名, Cr.pickleが出力される
                          cut_off=4.0, # cut_off
                          shuffle=True, # シャッフルするか
                          seed=1, # seed値
                          test_size=0.2, # test用の割合
                          test_output_dir="/nfshome17/knakajima/work/dataset/test/", # test用が出力されるパス
                          exclude_too_small_cell=True, # Trueならばcut_offの2倍以下のcellを持つ構造を除外する
                          exclude_too_large_force=True, # Trueならばforceが大きすぎる構造を除外する
                          max_allowable_force=50.0, # これ以上大きいforceを含む構造を除外する
                          ) 
```

## export_lammps_dumpposes
sfs.sfの構造を一つのfileにまとめる.(lammps形式のdumppos)
```python3
sfs.export_lammps_dumpposes(ofn="md.pos",  # 出力されるファイルのパス
                            out_columns=['type', 'x', 'y', 'z']) # 出力されるcolumns
```
<a id="anchor6"></a>
# analyze_frames

## count_mols
frameごとの分子数を数える.<br>
数えた結果はpandas.DataFrameとして得られる.
```python3
df_count_mols = sfs.count_mols(mode="bond_length", bond_length=[[1.2, 2.0],[2.0, 2.3]])
```
## count_bonds
frameごとの結合数を数える.
数えた結果はpandas.DataFrameとして得られる.
```python3
df_count_bonds = sfs.count_bnods(mode="bond_length", bond_length=[[1.2, 2.0],[2.0, 2.3]])
```
