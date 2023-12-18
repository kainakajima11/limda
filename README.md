# LiMDA
Library of Molecular Dynamics Analysis<br><br>
分子動力学シミュレーションにおける、ある一瞬の情報が入ったSimulationFrame(sf),<br>
sfを複数連ねたSimulationFrames(sfs)を用いた解析ライブラリ
## What you can do with LiMDA
- 分子動力学シミュレーションに関するファイルのimport, export
- VASP, Laich, Packmolといった分子動力学シミュレーションに必要な計算の実行
- 分子動力学シミュレーションの結果の解析
- ニューラルネットワークポテンシャルの学習のためのデータセットの作成
## Requirements
```
numpy>=1.24.3
pandas>=2.0.1
tqdm>=4.65.0
Cython=0.29.24
```
## Install
ソースコードをコピーした後、Cythonコードを用いるためのコンパイルが必要です。
```
git clone git@github.com:kainakajima11/limda.git
cd limda/src/limda
python3 setup.py build_ext --inplace
```
## Add Path
PYTHON PATH に $HOME/limda/src を追加してください。
```
#bashの場合 -> $HOME/.bashrcに下記を追記
export PYTHONPATH="$PYTHONPATH:$HOME/limda/src"
export PATH="${PATH}:$HOME/limda/src/scripts" # to use limda scripts
```
## How to update
各自ブランチを作りpull requestをしてください。<br>
その際、docstringを書いてください。

## Setting of VSCode
補完が出るようになります。<br>
nfshome17/knakajima を各自のものに変えてください。
```
"python.analysis.extraPaths": [
    "/nfshome17/knakajima/limda/src"
],
```
## Documents of limda
[SimulationFrame](https://github.com/kainakajima11/limda/blob/main/docs/SimulationFrame.md) <br>
[SimulationFrames](https://github.com/kainakajima11/limda/blob/main/docs/SimulationFrames.md) <br>
