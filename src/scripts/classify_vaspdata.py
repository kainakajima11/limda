#! /usr/bin/env python3

from pathlib import Path
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints
from pymatgen.io.vasp.outputs import Outcar, Oszicar, Vasprun
import os
from datetime import datetime
import shutil
from typing import Union
import argparse

class vaspDataManager:
    classification: list[str]
    info: list[str]
    warning: bool
    warn_sentence: list[str]

    def __init__(self,):
        self.info = {}
        self.classification = []
        self.warning = False
        self.warn_sentence = []

    def check_poscar(self, path):
        data = Poscar.from_file(path)
        self.classification.append(f"{''.join(data.site_symbols)}_")
        self.info["ATOM_INFO"] = data.structure.composition.formula
        cell = data.structure.lattice
        self.info["CELL"] = f"{cell.a} {cell.b} {cell.c}"

    def check_incar(self, path):
        data = Incar.from_file(path).as_dict()
        # Ensumble
        if "ISIF" in data:
            if data["ISIF"] == 3:
                self.classification.append("NPT_")
                if "PSTRESS" in data:
                    self.classification.append(f"{int(data['PSTRESS'])}kbar_")                    
        # Temperature
        if "TEBEG" in data and "TEEND" in data:
            if int(data["TEBEG"]) == int(data["TEEND"]):
                self.classification.append(f"{int(data['TEBEG'])}K_")
            else:
                self.classification.append(f"{int(data['TEBEG'])}_{int(data['TEEND'])}K_")       
        # SMEARING
        if "ISMEAR" in data and data["ISMEAR"] == 1:
            self.classification.append("SMR_")
        # IVDW DFT-D3
        if "IVDW" in data and data["IVDW"] != 0:
            self.classification.append(f"D3_{data['IVDW']}_")
        # DFT+U 
        if "LDAU" in data and data["LDAU"] == ".TRUE.":
            self.classification.append("+U_")    
        # SPIN
        if "ISPIN" in data and data["ISPIN"] == 2:
            self.classification.append("SPIN_")
        # ENCUT
        if "ENCUT" in data:
            self.classification.append(f"ECT_{data['ENCUT']}_")             
        # EDIFF
        if "EDIFF" in data and data["EDIFF"] > 1e-5:
            self.warning = True
            self.warn_sentence.append(f"EDIFF IS LARGER THAN 1e-5, {data['EDIFF']}") 
        # PREC
        if "PREC" not in data:
            self.warning = True
            self.warn_sentence.append("PREC IS NOT Accurate") 
        elif data["PREC"] != "Accurate":
            self.warning = True
            self.warn_sentence.append(f"PREC IS NOT Accurate {data['PREC']}") 

    def check_kpoints(self, path):
        data = Kpoints.from_file(path)
        if any(kpt != 1 for kpt in data.kpts[0]):
            self.classification.append(f"KPT_{''.join(data.kpts[0])}")

    def check_oszicar(self, path):
        data = Oszicar(path)
        self.classification.append(f"{len(data.ionic_steps)}stps_")           

    def check_files(self, path):
        self.check_poscar(path / "POSCAR")
        self.check_incar(path / "INCAR")
        self.check_kpoints(path / "KPOINTS")
        self.check_oszicar(path / "OSZICAR")
        if self.warning:
            self.classification.insert(0, "Warning_")
        
class Input:
    from_dir_path: Union[str, Path]
    vaspdata_paths: list[str]
    to_dir_name: str

    def __init__(self, from_dir_path, to_dir_name,):
        self.from_dir_path = from_dir_path
        self.to_dir_name = to_dir_name
        self.vaspdata_paths = []
    
    def search_vaspdirs(self):
        outcar_paths = list(Path(self.from_dir_path).glob(f"**/OUTCAR"))
        for path in outcar_paths:
            self.vaspdata_paths.append(path.parent)

    def copy_vaspfile(self, path, new_dir_path, files):
        for file in files:
            shutil.copyfile(Path(path) / file, Path(new_dir_path) / file)
        
def deal_target_dir(info):
    assert(len(info) == 3)
    ipt = Input(info[1], info[2])
    return ipt

def init_info(path) -> (str, list[Input]):
    storage_path = None
    inputs = []
    move_files = []
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "#target_dir":
                ipt = deal_target_dir(spline)
                ipt.search_vaspdirs()
                inputs.append(ipt)
            if spline[0] == "#storage_dir":
                storage_path = spline[1]
            if spline[0] == "#move_file":
                for i in range(1,len(spline)):
                    move_files.append(spline[i])
    return (storage_path, inputs, move_files)

def add_basic_info(path):
    with open(path, "a") as f:
        now_str = datetime.now().strftime("%Y-%m-%d-%H:%M")
        f.write(f"Following Data were Stored in {now_str} by {os.environ.get('USER')}\n\n")

def add_readme(readme_path, manager, new_dir_path, path):
    with open(readme_path, "a") as f:
        f.write(f"---- {new_dir_path.name} ----\n\n")
        if manager.warning:
            f.write("-#-#-#-# WARNING #-#-#-#-\n")  
            f.write("\n".join(manager.warn_sentence))
            f.write("\n-#-#-#-#-#-#-#-#-#-#-#-#-\n\n")
        f.write(f"COMPOSITION : {manager.info['ATOM_INFO']}\n")
        f.write(f"CELL LATTICE : {manager.info['CELL']}\n")
        f.write(f"from {path}\n\n\n")

def classify_and_move_vasp_dir(storage_path: Union[str, Path], ipt: Input, move_files: int):
    storage_path = Path(storage_path)
    dir_i = 0
    if not os.path.exists(storage_path / ipt.to_dir_name):
        os.mkdir(storage_path / ipt.to_dir_name)
    else:
        dir_i = len(os.listdir(storage_path / ipt.to_dir_name))

    add_basic_info(storage_path / ipt.to_dir_name / "README")    

    managers = [vaspDataManager() for _ in range(len(ipt.vaspdata_paths))]
    for i, path in enumerate(ipt.vaspdata_paths):
        managers[i].check_files(path)

        new_dir_path = storage_path / ipt.to_dir_name / str('{:03d}'.format(dir_i) + "_" + ''.join(managers[i].classification).rstrip("_"))
        os.mkdir(new_dir_path)

        ipt.copy_vaspfile(path, new_dir_path, move_files)

        add_readme(storage_path / ipt.to_dir_name / "README", managers[i], new_dir_path, path)

        dir_i += 1

        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="VASPのdataを分類・移動します")
    parser.add_argument("input_file_name", default=None, type=str, help="情報が入ったファイル")
    args = parser.parse_args()
    assert args.input_file_name is not None, "input_file_nameを設定してください"

    storage_path, inputs, move_files = init_info(args.input_file_name)
    for ipt in inputs:
        classify_and_move_vasp_dir(storage_path, ipt, move_files)