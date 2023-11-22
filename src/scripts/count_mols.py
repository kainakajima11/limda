#! /usr/bin/env python3
from limda import SimulationFrame, SimulationFrames
import pandas as pd
import argparse
from pprint import pprint

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="分子数をカウントする")
    parser.add_argument(
        "-i",
        "--input_file_name",
        default=None,
        type=str,
        help="対象のファイル名、指定しない場合はカレントディレクトリのdump.posファイルに対して分子をカウントする",
    )
    parser.add_argument(
        "-p", "--para_str", default="", type=str, help="para, ex. 'Cr Mn Fe Co Ni'"
    )
    parser.add_argument(
        "-s", "--skip_num", default=1, type=int, help="何個おきにdumpposを読み込むか"
    )

    args = parser.parse_args()

    if args.input_file_name is not None:
        sf = SimulationFrame()
        sf.import_para_from_str(args.para_str)
        sf.import_file(args.input_file_name)
        pprint(sf.count_mols())
    else:
        sfs = SimulationFrames()
        sfs.import_dumpposes(skip_num=args.skip_num)
        df = sfs.count_mols()
        print(df)
