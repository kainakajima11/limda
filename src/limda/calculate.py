import numpy as np
import pathlib
from tqdm import tqdm
import subprocess
import os
import time
from typing import Union

try:
    from ase.build import molecule
    from ase import Atoms
except:
    pass

class Calculate(

):
    """計算を実行するクラス
    """
    def __init__(self):
        pass

    def vasp(
            self,
            calc_directory: str,
            poscar_comment: str="",
            poscar_scaling_factor: float=1.0,
            incar_config: dict=None,
            potcar_root: str=None,
            kpoints_comment: str="",
            kpoints_kx: int=1,
            kpoints_ky: int=1,
            kpoints_kz: int=1,
            iconst_config: list[str]=None,
            vasp_command: str="vasp_std",
            print_vasp: bool=True,
            exist_ok: bool=False,
            poscar_from_contcar: bool = False,
            contcar_path: str = ""
    ):
        """vaspを実行する.
        Parameters
        ----------
            calc_directory: str
                vaspで計算を実行するディレクトリ
            poscar_comment: str
                POSCARの1行目に書かれるコメント
            poscar_scaling_factor: float
                vaspのドキュメントをを参照してください.基本的に変える必要はない
            potcar_root: str
                元のPOTCARがあるフォルダ
            incar_config: dict
                INCARに書かれるconfig
            kpoints_comment: str
                KPOINTSの1行目に書かれるコメント
            kpoints_kx: int
            kpoints_ky: int
            kpoints_kz: int
                x, y, z方向のK点     
            iconst_config: List[str]
                ICONSTに書かれるconfig
                VASPでNPT計算をする時で、セルの角度を固定したい時は以下のように指定する
                iconst_config = ['LA 1 2 0',
                                'LA 1 3 0',
                                'LA 2 3 0']
            vasp_command: str
                実行するvaspのコマンド
            print_vasp: bool
                実行中のvaspの標準出力をpythonからも表示するかどうか
            exist_ok: bool
                calc_directoryが存在する時、exist_ok=Trueのときは上書きしてvaspを実行する
                calc_directoryが存在する時、exist_ok=Falseのときは上書きしない
            poscar_from_contcar: bool
                poscarをcontcarから作るときTrueにする
            contcar_path: str
                contcarからposcarを作るときに必要なcontcarのpath
        """
        assert potcar_root is not None, "There isn't potcar_root."
        assert incar_config is not None, "There isn't incar_config."
        calc_directory = pathlib.Path(calc_directory)
        for key in ["NCORE", "NPAR", "KPAR"]:
            if key not in incar_config:
                incar_config[key] = 1

        os.makedirs(calc_directory, exist_ok=exist_ok)
        poscar_path = calc_directory / "POSCAR"
        incar_path = calc_directory / "INCAR"
        potcar_path = calc_directory / "POTCAR"
        kpoints_path = calc_directory / "KPOINTS"
        iconst_path = calc_directory / "ICONST"
        num_process = incar_config["NCORE"] * incar_config["NPAR"] * incar_config["KPAR"]
        if not poscar_from_contcar:
            self.export_vasp_poscar(poscar_path, poscar_comment, poscar_scaling_factor)
        else:
            self.export_vasp_poscar_from_contcar(poscar_path, contcar_path)
        self.export_vasp_incar(incar_path, incar_config)
        self.export_vasp_potcar(potcar_path, potcar_root)
        self.export_vasp_kpoints(kpoints_path, kpoints_comment, kpoints_kx, kpoints_ky, kpoints_kz)

        if iconst_config is not None:
            self.export_vasp_iconst(iconst_path, iconst_config)

        cmd = f'mpiexec -np {num_process} {vasp_command} > stdout'
        vasp_md_process = subprocess.Popen(cmd, cwd=calc_directory, shell=True)
        time.sleep(5)
        if print_vasp:
            tail_process = subprocess.Popen(f'tail -F stdout', cwd=calc_directory, shell=True)
            while vasp_md_process.poll() is None:
                time.sleep(1)
            tail_process.kill()
#----------------------------------------------------------------------------------------------
    def laich(self,
              calc_dir: str='laich_calc',
              para_file_path: str=None,
              laich_cmd: str ='laich',
              laich_config :dict=None,
              print_laich: bool=False,
              exist_ok=False):
        """LaichでMD,または構造最適化を実行する。
        Parameters
        ----------
            calc_dir: str
                laichの結果が出力される部分
            para_file_path: str
                para fileのpath
            laich_mode: str
                "MD" or "OPT" にする.
                MD : 分子動力学計算を行う.
                OPT : 構造最適化を行う.
            laich_config: dict
                変数の入ったdict
            print_laich: bool
                out fileをprintするか
            exist_ok: bool
        """
        calc_dir = pathlib.Path(calc_dir)
        if not exist_ok:
            assert not calc_dir.exists(), "calc_dir exists."
        calc_dir.mkdir()

        if para_file_path is not None:
            para_file_path = pathlib.Path(para_file_path)
        input_file_path = calc_dir / 'input.rd'
        config_file_path = calc_dir / 'config.rd'
        self.export_input(input_file_path)
        with open(config_file_path, 'w') as f:
            for key in laich_config.keys():
                f.write(f"{key} {laich_config[key]}\n")

        if para_file_path is not None:
            try:
                subprocess.run(f'cp {para_file_path} {calc_dir / "para.rd"}', shell=True)
            except:
                pass

        np = laich_config["MPIGridX"]*laich_config["MPIGridY"]*laich_config["MPIGridZ"]

        cmd = f"mpiexec.hydra -np {np} {laich_cmd} < /dev/null >& out"
        laich_process = subprocess.Popen(cmd, cwd=calc_dir, shell=True)
        time.sleep(5)
        if print_laich:
            tail_process = subprocess.Popen(f'tail -F out', cwd=calc_dir, shell=True)
            while laich_process.poll() is None:
                time.sleep(1)
            tail_process.kill()

        dumppos_paths = list(calc_dir.glob('./dump.pos.*'))
        dumppos_paths.sort(reverse=True)
        assert len(dumppos_paths) != 0, "dumpposが生成されていません"
        optimized_dumppos_path = dumppos_paths[0]
        self.import_dumppos(optimized_dumppos_path)
#------------------------------------------------------------------------------------------
# laich_config = {
#     "Mode": "MD"
#     "ForceField": "Reaxff",
#     "XYZFile": "input.rd",
#     "ParaFile": "para.rd",
#     "TimeStep": 0.25,
#     "TotalStep": 10000,
#     "ObserveStep": 1,
#     "FileStep": 1000,
#     "BondStep": 1000,
#     "SaveRestartStep": 10000,
#     "NPUGS": 1,
#     "NNPModelPath": "script_model.pth",
#     "WEIGHTPATH": './script_model.pth',
#     "MPIGridX": 1,
#     "MPIGridY": 1,
#     "MPIGridZ": 1,
#     "CUTOFF": 10.0,
#     "MARGIN": 1.0,
#     "GhostFactor": 20.0,
#     # MD Mode
#     "NNPModelPath": "script_model.pth",
#     "OMPGridX": 1,
#     "OMPGridY": 1,
#     "OMPGridZ": 1,
#     "ShowMask": 1,
#     "ReadVelocity": 0,
#     "Thermo": "Langevin",
#     "AimTemp": 300.0,
#     "InitTemp": 300.0,
#     "FinalTemp": 300.0,
#     "ThermoFreq": 0.005,
#     # OPT Mode
#     "DelR": 0.0001,
#     "MaxR": 0.1
# }
#---------------------------------------------------------------------------------------------
    def packmol(self,
                sf_list: list,
                pack_num_list: list[int],
                tolerance: float=2.0,
                packmol_tmp_dir: Union[str,pathlib.Path]="./packmol_tmp",
                xyz_condition: list[float]=None,
                seed: int=-1,
                print_packmol=False
                ):
        """packmolで原子を詰める
        Parameters
        ----------
            sf_list : list[SimulationFrame]
                詰めるsfのlist
            pack_num_list : list
                詰めるsdatの個数のリスト
            tolerance : float
                最小の原子間距離, 原子を詰める時に原子間がtolerance以上になるように詰める
            xyz_condition : list
                詰める原子の座標の条件
            packmol_tmp_dir : Union[str,pathlib.Path]
                packmolを動かすときの一時的なディレクトリ
            seed : int
                シード値
                seed = -1のときはseedは時間で決定される
            print_packmol : bool
                packmolの標準出力を表示するか
        Example
        -------
        xyz_condition = [
          # [xmin, ymin, zmin, xmax, ymax, zmax] で指定する
            [2, 2, 2, 8, 18, 28],  # sf2の条件
            [2, 2, 2, 8, 10, 10]   # sf3の条件
        ]
        sf1.packmol(sf_list=[sdat2, sdat3], pack_num_list=[5, 8], xyz_condition=xyz_condition)
        とすると, sf1にsf2の原子が5個, sf3の原子が8個詰められる
        sf2は2 <= x <= 8 かつ 2 <= y <= 18 かつ 2 <= z <= 28 の位置のみに詰められる
        sdat3は2 <= x <= 8 かつ 2 <= y <= 10 かつ 2 <= z <= 10 の位置のみに詰められる

        Note
        ----
            周期境界条件の場合は端まで詰めると計算が回らなくなるので注意.
            境界には間を開けることを推奨

            この関数を使うには、$ packmol のみでコマンドが使えるようにパスを通しておく必要がある
        """
        if xyz_condition is not None:
            assert len(xyz_condition) == len(sf_list), "sf_listとxyz_conditonを同じ長さにしてください."
        
        packmol_tmp_dir = pathlib.Path(packmol_tmp_dir)
        packmol_tmp_dir.mkdir()

        first_line = [
            f"tolerance {tolerance}\n",
            f"filetype xyz\n",
            f"seed {seed}\n",
            f"output packmol_mixture_result.xyz\n\n"
        ]

        second_line = []
        if self.atoms is not None and self.get_total_atoms() >= 1:
            second_line = [
                "structure this_sf.xyz\n",
                "\tnumber 1\n",
                "\tfixed 0. 0. 0. 0. 0. 0. 0. \n",
                "end structure\n\n",
            ]

        with open(packmol_tmp_dir / "packmol_mixture_comment.inp", 'w') as f:
            f.writelines(first_line)
            f.writelines(second_line)
            for frame_idx in range(len(sf_list)):
                f.write(f"structure sf_idx_{frame_idx}.xyz\n")
                f.write(f"\tnumber {pack_num_list[frame_idx]}\n")
                if xyz_condition is not None:
                    f.write(f"\tinside box")
                    for condition in xyz_condition[frame_idx]:
                        f.write(f"\t{condition}")
                    f.write("\n")
                f.write("end structure\n\n")
        
        # export xyz file
        if self.atoms is not None and self.get_total_atoms() >= 1:
            self.export_xyz(packmol_tmp_dir / "this_sf.xyz", structure_name="this_sf")

        for frame_idx in range(len(sf_list)):
            sf_list[frame_idx].export_xyz(packmol_tmp_dir / f"sf_idx_{frame_idx}.xyz",
                                          structure_name=f"sf_idx_{frame_idx}")
        
        # run packmol
        cmd = f"packmol < {'packmol_mixture_comment.inp'}"
        if print_packmol:
            p = subprocess.Popen(cmd, cwd=packmol_tmp_dir, shell=True, )
        else:
            p = subprocess.Popen(cmd, cwd=packmol_tmp_dir, shell=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        p.wait()
        # import result
        self.import_xyz(packmol_tmp_dir / "packmol_mixture_result.xyz")


