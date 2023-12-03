#! /usr/bin/env python3
from limda import SimulationFrames
import pandas as pd
import argparse

pd.set_option("display.max_rows", None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="dumppos fileをひとまとめにする")
    parser.add_argument(
        "dir_path", default=None, type=str, help="dumppos fileがあるdirectory"
    )
    parser.add_argument(
        "-o", "--output_file_name", default="md.pos", type=str, help="dumpposをまとめたpos fileの名前"
    )
    parser.add_argument(
        "-p", "--para_str", default="", type=str, help="para, ex. 'Cr Mn Fe Co Ni'"
    )
    parser.add_argument(
        "-s", "--skip_num", default=None, type=int,  help="いくつおきにdumpposを読み込むのか"
    )
    args = parser.parse_args()

    assert args.dir_path is not None, "dir_pathを設定してください"

    sfs = SimulationFrames()
    sfs.import_para_from_str(args.para_str)
    sfs.import_dumpposes(dir_name=args.dir_path, skip_num=args.skip_num)
    sfs.export_lammps_dumpposes(args.output_file_name)