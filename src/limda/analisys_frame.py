from .neighbor import cy_get_neighbor_list
from .neighbor import cy_count_molecules

class AnalisysFrame(

):
    def __init__(self):
        pass

    def get_neighbor_list(self, cut_off: float=3.4 ,bond_length: list[list[float]] = [])->list[list[int]]:
        """neighbor_list を作成します。
        Parameters
        -----------
            cut_off: float
                typeに関係なく、距離を指定するときに用いる。
                bond_lengthを指定しなければ、この値が採用されます。
            bond_length: list[list[float]]
                size : type数 x type数
                Example
                -------
                "C H O" というparaならば bond_length[0][1] : C-H の最大結合距離 = bond_length[1][0]
        Return val
        ----------
            neighbor_list: list[list[int]]
                list[i番目の原子と結合している原子のidが入ったlist]
        """
        if not bond_length:
            if cut_off*3 > min(self.cell):
                mesh_length = min(self.cell)/3-0.01
            else:
                mesh_length = cut_off
            neighbor_list = cy_get_neighbor_list(atoms_type = self.atoms["type"],
                                                 atoms_pos = [self.atoms["x"], self.atoms["y"], self.atoms["z"]],
                                                 mesh_length = mesh_length + 0.01,
                                                 atom_num = len(self),
                                                 bond_length = [],
                                                 cut_off = cut_off,
                                                 cell = self.cell,
                                                 mode = "neighbor")
        else:
            mesh_length = max(list(map(lambda x: max(x), bond_length)))
            if mesh_length*3 > min(self.cell):
                mesh_length = min(self.cell)/3-0.01
            neighbor_list = cy_get_neighbor_list(atoms_type = self.atoms["type"],
                                                atoms_pos = [self.atoms["x"], self.atoms["y"], self.atoms["z"]],
                                                mesh_length = mesh_length + 0.01,
                                                atom_num = len(self),
                                                bond_length = bond_length,
                                                cut_off = 0,
                                                cell = self.cell,
                                                mode = "neighbor")
        for idx in range(len(neighbor_list)):
            neighbor_list[idx] = sorted(neighbor_list[idx])

        return neighbor_list
#-------------------------------------------------------------------------------------
    def get_edge_idx(self, cut_off: float)->list[list[int]]:
        """allegroのデータセットを作るように,edge_indexを作成します。
            neighbor_listとは、listのサイズが異なり、list[list[int]] : 2 x 結合個数で重複はなしです。
            Parameters
            ----------
            cut_off: float
                edgeとしてみなす最大距離
        """
        if cut_off*3 > min(self.cell):
            mesh_length = min(self.cell)/3-0.01
        else:
            mesh_length = cut_off
        edge_idx = cy_get_neighbor_list(atoms_type = self.atoms["type"],
                                   atoms_pos = [self.atoms["x"],self.atoms["y"], self.atoms["z"]],
                                   mesh_length = mesh_length+0.01,
                                   atom_num = len(self),
                                   bond_length = [],
                                   cut_off = cut_off,
                                   cell = self.cell,
                                   mode = "edge")
        return edge_idx
        
#-------------------------------------------------------------------------------------    
    def get_neighbor_list_test(self, bond_length: list[list[float]])->list[list[int]]:
        """ neighbor_listを作成します。
            pythonでO(N^2)のため、get_neighbor_list()のtest用です。
        Parameters
        -----------
            bond_length: list[list[float]]
        """
        neighbor_list_test = [[] for _ in range(len(self))]
        x = self.atoms['x'].values
        y = self.atoms['y'].values
        z = self.atoms['z'].values
        atom_types = self.atoms['type'].to_list()
        for i  in range(self.get_total_atoms()):
            for j in range(i+1, self.get_total_atoms()):
                dx:list[float] = [None,None,None]
                dx[0] = x[j] - x[i]
                dx[1] = y[j] - y[i]
                dx[2] = z[j] - z[i]
                for ax in range(3):
                    if dx[ax] < -self.cell[ax]/2:
                        dx[ax] += self.cell[ax]
                    elif self.cell[ax]/2 < dx[ax]:
                        dx[ax] -= self.cell[ax]
                if dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] <= bond_length[atom_types[i]-1][atom_types[j]-1]*bond_length[atom_types[i]-1][atom_types[j]-1]:
                    neighbor_list_test[i].append(j)
                    neighbor_list_test[j].append(i)

        for idx in range(len(neighbor_list_test)):
            neighbor_list_test[idx] = sorted(neighbor_list_test[idx])

        return neighbor_list_test
    
#--------------------------------------------------------------------------------------    
    def count_molecules(self, cut_off: float=3.4 ,bond_length: list[list[float]] = []):
        """
        系内に何の分子が何個存在するか数え上げます。
        bond_lengthを引数で指定した場合、bond_lengthをもとに数えます。
        指定しなければcut_offをもとに数えます。
        Parameters
        ----------
        cut_off: float
            結合距離をすべて同じで探索するときに指定する結合距離
        bond_length: list[list[float]]
            原子タイプごとに結合距離を指定したいときに使用します。
        
        Example
        -------
            系にCO2が3つ、H2Oが1つ存在するとき
            
            C1O2 : 3
            H2O1 : 1

            とprintされます。
        """
        if not bond_length:
            if cut_off*3 > min(self.cell):
                mesh_length = min(self.cell)/3-0.01
            else:
                mesh_length = cut_off
            molecules_list = cy_count_molecules(atoms_type = self.atoms["type"],
                                                 atoms_pos = [self.atoms["x"], self.atoms["y"], self.atoms["z"]],
                                                 mesh_length = mesh_length + 0.01,
                                                 atom_num = len(self),
                                                 bond_length = [],
                                                 cut_off = cut_off,
                                                 cell = self.cell,
                                                 typ_len = len(self.atom_symbol_to_type))
        else:
            mesh_length = max(list(map(lambda x: max(x), bond_length)))
            if mesh_length*3 > min(self.cell):
                mesh_length = min(self.cell)/3-0.01
            molecules_list = cy_count_molecules(atoms_type = self.atoms["type"],
                                                atoms_pos = [self.atoms["x"], self.atoms["y"], self.atoms["z"]],
                                                mesh_length = mesh_length + 0.01,
                                                atom_num = len(self),
                                                bond_length = bond_length,
                                                cut_off = 0,
                                                cell = self.cell,
                                                typ_len = len(self.atom_symbol_to_type))
            
        molecules_dict :dict[str, int] = {}
        for molecule in molecules_list:
            mol_str = ""
            for i,atom_number in enumerate(molecule):
                if not atom_number:
                    continue
                if mol_str:
                    mol_str = f"{mol_str} "
                mol_str = f"{mol_str}{self.atom_type_to_symbol[i+1]}{atom_number}"
            if mol_str not in molecules_dict.keys():
                molecules_dict[mol_str] = 1
            else:
                molecules_dict[mol_str] += 1

        for molecule, num in molecules_dict.items():
            print(f"{molecule} : {num}")