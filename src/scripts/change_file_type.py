git #! /usr/bin/env python3
from limda import SimulationFrame
import pandas as pd
import argparse

pd.set_option('display.max_rows', None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="ファイルの変換"
    )
    parser.add_argument("-f", "--file_name", default=None, type=str,
                        help = "変換前のファイル名")
    parser.add_argument("-o", "--output_file_name", default=None, type=str,
                        help = "変換後のファイル名")    
    parser.add_argument("-p", "--para_str", default=None, type=str,
                        help = "para, ex. 'Cr Mn Fe Co Ni'")
    
    args = parser.parse_args()

    assert args.para_str is not None,"para_strを設定してください"
    assert args.file_name is not None,"file_nameを設定してください"
    assert args.output_file_name is not None,"output_file_nameを設定してください"

    sf = SimulationFrame()
    sf.import_para_from_str(args.para_str)
    sf.import_file(args.file_name)
    sf.export_file(args.output_file_name)
    