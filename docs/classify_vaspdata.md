# classify_vaspdata.py
ターミナル上で実行できるスクリプトです。<br>
VASPを回した後のディレクトリに対して、分類をしたり、別の場所に保管するために使います。<br>
POSCAR, INCAR, OSZICAR, KPOINTS から情報を抜き出します。<br>

### 実行の方法
```
classify_vaspdata.py input
```
第一引数にinput fileを指定する必要があります。<br>

### input fileの書き方
#### 例
```
#storage_dir /nfshome17/vaspdata_storage
#move_file INCAR OUTCAR POSCAR POTCAR OSZICAR CONTCAR
#target_dir /nfshome17/knakajima/work/fpc/B/FCC/B_F_L Cantor_FCC_bulk
#target_dir /nfshome17/knakajima/work/fpc/OH/H2O_density1.0 H2O
```
#### #target_dir
#target_dir 第一引数 第二引数 <br>
第一引数 : 対象となるディレクトリのパス
第二引数 : storage_dirにできる新しいディレクトリの名前<br><br>
複数のディレクトリに対して操作がしたければ、#target_dirの行を増やしてください<br>
第二引数にはなるべく多くの情報をわかりやすくした名前を入れてください。<br><br>
例
```
Cantor+H2O_FCC_slab_013 
```
特に構造に関しては、分類が難しいためどんな構造かわかるようにしてください。<br>
##### 基本的に含ませる情報
 - 何の材料か (Cantor+H2O, MoDTC+Fe ...)
 - 結晶の構造 (FCC, BCC, Amorphous)
 - 表面の情報 (slab_001, bulk)

##### 他の情報の例
 - 初期速度を含む (with_velocity)
 - 不安定な構造 (unstable)
 - 格子定数が特殊 (lattice_const_3.5, strain)

#### #storage_dir 
#storage_dir 第一引数<br>
第一引数 : 指定したファイルが保管されるパス<br>

#### #move_file
#move_file 第一引数 第二引数 ... <br>
第n引数 : 保管するファイル名<br>
ディレクトリの中の何のファイルを保管するか<br>

### 使用例
対象となるvasp dataのディレクトリ構造
```
#target_dir 1 argument
├── 〇〇_0
│   ├── #move_file 1 argument
│   ├── #move_file 2 argument
│   ├── #move_file 2 argument
│   ├── #move_file 3 argument
│   ├── file
│   └── file
├── 〇〇_1
│   ├── #move_file 1 argument
│   ├── #move_file 2 argument
│   ├── #move_file 2 argument
│   ├── #move_file 3 argument
│   ├── file
│   └── file 
.   .
.   .
.   .
``` 

実行後のstorage_dirのディレクトリ構造<br>
VASP実行ディレクトリがそれぞれ分類されます。(000..., 001..., )<br>
```
#storage_dir 1 argument
├── #target_dir 2 argument
│   ├── 000...
|   │   ├── #move_file 1 argument
|   │   ├── #move_file 2 argument
|   │   ├── #move_file 2 argument
|   │   └── #move_file 3 argument
│   ├── 001...
|   │   ├── #move_file 1 argument
|   │   ├── #move_file 2 argument
|   │   ├── #move_file 2 argument
|   │   └── #move_file 3 argument
│   ├── 002...
|   │   ├── #move_file 1 argument
|   │   ├── #move_file 2 argument
|   │   ├── #move_file 2 argument
|   │   └── #move_file 3 argument
.   .    
.   .   
.   └── README
```
### 分類方法
ディレクトリの名前に追加される、もしくはREADMEに書かれます。<br>
```
003_CrMnFeCoNi_NPT_10kbar_600K_SMR_SPIN_ECT_400_5stps
```
##### 分類の種類(ディレクトリの名前になる)
 - インデックス (000, 001, ..., 999)
 - 構成元素 (CrMnFeCoNi, HO)
 - アンサンブル NVT以外なら表記(NPT_10kbar)
 - 温度 (600K, 300-2300K)
 - スメアリング 有効なら表記(SMR)
 - スピン 有効なら表記(SPIN)
 - ENCUT 設定していれば表記(ECT_400)
 - step数 (5stps, 400stps)
 - DFT+U 有効なら表記(+U)
 - DFT-D3 有効なら表記(D3_12)
 - K点 111以外なら表記(KPT_222)

### README
細かい情報や分類を行った時の状況をREADMEに記述します。<br>
例<br>
```
Following Data were Stored in 2024-04-08-11:30 by knakajima

---- 000_MnFeCoNi_NPT_10kbar_600K_SMR_SPIN_ECT_400_5stps ----

COMPOSITION : Mn43 Fe21 Co21 Ni23
from /nfshome17/knakajima/work/fpc/B/FCC/B_F_L/B_F_L_111
```
 - 実行日時
 - 実行した人
 - 組成 (ディレクトリごと)
 - どこからファイルを持ってきたか (ディレクトリごと）

