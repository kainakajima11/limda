#! /usr/bin/env python3
from limda import SimulationFrame
import pandas as pd
import argparse

pd.set_option("display.max_rows", None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="密度を求める")
    parser.add_argument(
        "-i", "--input_file_name", default=None, type=str, help="対象のファイル名"
    )
    parser.add_argument(
        "-p", "--para_str", default="", type=str, help="para, ex. 'Cr Mn Fe Co Ni'"
    )

    args = parser.parse_args()

    assert args.input_file_name is not None, "file_nameを設定してください"

    sf = SimulationFrame()
    sf.import_para_from_str(args.para_str)
    sf.import_file(args.input_file_name)
    print(sf.density())
