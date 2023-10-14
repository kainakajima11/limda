#! /usr/bin/env python3

import argparse
import pathlib
import yaml
from natsort import natsorted
from limda import SimulationFrames

"""outcar_to_dumppos.py
vaspのOUTCARファイルをlammpsのdumpposファイルに変換するpythonスクリプト
OUTCARファイルがあるフォルダにmd.posを作成する
```sh
outcar_to_dumppos.py path_to_outcar_dir para_list
```
で実行できる
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='vaspのOUTCARファイルをlammpsのdumpposファイルに変換する'
    )
    parser.add_argument("target_dir", help="OUTCARファイルがあるフォルダ")
    parser.add_argument('para_list', help='原子のリスト.例えば、C H O N Si など', nargs='+', type=str)
    args = parser.parse_args()

    target_dir = pathlib.Path(args.target_dir)
    para_list = args.para_list

    print(f"directory : {str(target_dir.resolve())}")
    print(f"para : {para_list}")

    outcar_paths = natsorted(list(target_dir.glob("**/OUTCAR")))
    sfs_list = []

    for outcar_path in outcar_paths:
        vasp_dir = outcar_path.parent
        print(vasp_dir)

        sfs = SimulationFrames()
        sfs.import_para_from_list(para_list)
        sfs.import_vasp(vasp_dir)
        if not len(outcar_paths) == 1:
            sfs.export_lammps_dumpposes(vasp_dir / "md.pos")

        sfs_list.append(sfs)

    sfs_output_all_frames = SimulationFrames()
    sfs_output_all_frames.import_para_from_list(para_list)
    sfs_output_all_frames.concat_sfs(sfs_list)
    sfs_output_all_frames.export_lammps_dumpposes(target_dir / "md.pos")
