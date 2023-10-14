#! /usr/bin/env python3

import argparse
import pathlib
import yaml
from limda import SimulationFrames

"""make_allegro_dataset.py
vaspのOUTCARファイルをallegro_datasetに変換するスクリプト
```sh
make_allegro_dataset.py dataset_config.yaml
```
で実行できる
"""

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description="vasp:OUTCARをallegroのdatasetに変換"
    )
    parser.add_argument("dataset_config", help="dataset作成のためのconfig")

    args = parser.parse_args()

    dataset_config = pathlib.Path(args.dataset_config)

    with open(dataset_config, "r") as f:
        config = yaml.safe_load(f)

    assert "cut_off" in config
    assert "dataset_outcar_dirs" in config
    assert len(config["dataset_outcar_dirs"]) >= 1
    assert "para_list" in config
    assert "train_dataset_dir" in config
    assert "shuffle_frames" in config
    assert "test_size" in config
    if config["test_size"] != 0.0:
        assert "test_dataset_dir" in config

    config["train_dataset_dir"] = pathlib.Path(config["train_dataset_dir"])
    if config["test_size" ] != 0.0:
        config["test_dataset_dir"] = pathlib.Path(config["test_dataset_dir"])

    config["dataset_outcar_dirs"] = [pathlib.Path(directory) for directory in config["dataset_outcar_dirs"]]
    for outcar_dir in config["dataset_outcar_dirs"]:
        print(outcar_dir.resolve())
        outcar_paths = outcar_dir.glob("**/OUTCAR")

        for outcar_path in outcar_paths:
            vasp_dir = outcar_path.parent
            print(vasp_dir)

            sfs = SimulationFrames()
            sfs.import_para_from_list(config["para_list"])
            sfs.import_vasp(vasp_dir)
            
            if config['test_size'] == 0:
                sfs.export_allegro_frames(config['train_dataset_dir'] / vasp_dir.parent.name, 
                                        f"{vasp_dir.name}", 
                                        cut_off=config['cut_off'], 
                                        shuffle=config['shuffle_frames'], 
                                        )
            else:
                sfs.export_allegro_frames(config['train_dataset_dir'] / vasp_dir.parent.name, 
                                        f"{vasp_dir.name}", 
                                        cut_off=config['cut_off'], 
                                        shuffle=config['shuffle_frames'], 
                                        test_size=config['test_size'], 
                                        test_output_dir=config['test_dataset_dir'] / vasp_dir.parent.name, test_output_file_name=f"{vasp_dir.name}"
                                        )