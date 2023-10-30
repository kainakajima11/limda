import pandas as pd


class AnalyzeFrames:
    def __init__(self):
        pass

    def count_mols(
        self,
        mode: str = "bond_length",
        cut_off: float = None,
        bond_length: list[list[float]] = None,
    ) -> dict[str, int]:
        """分子数を数える
        例えば、水分子が3個とアンモニアが1個あるときは
        {"H2O1":3,
         "H3N1":1}
        Parameters
        ----------
            mode: str
                "bond_length"または"cut_off"
                mode = "bond_length"とした場合はneighbor listを結合種の長さ(bond_length)によって作成する
                mode = "cut_off"とした場合はneighbor listをカットオフによって作成する
            cut_off: float
                カットオフ半径
            bond_length: list[list[float]]
                結合の長さ
        """
        count_mols_lists = []
        for frame_idx in range(len(self.sf)):
            count_mols_lists.append(
                self.sf[frame_idx].count_mols(
                    mode=mode, cut_off=cut_off, bond_length=bond_length
                )
            )
        df_count_mols = pd.DataFrame(count_mols_lists).fillna(0).astype(int)
        df_count_mols.index = self.get_step_nums()

        return df_count_mols

    def count_bonds(
        self,
        mode: str = "bond_length",
        cut_off: float = None,
        bond_length: list[list[float]] = None,
    ) -> dict[str, int]:
        """結合数を数える
        例えば、水分子が3個あるときは
        {"H-O": 9, "H-H": 0, "O-O": 0}
        Parameters
        ----------
            mode: str
                "bond_length"または"cut_off"
                mode = "bond_length"とした場合はneighbor listを結合種の長さ(bond_length)によって作成する
                mode = "cut_off"とした場合はneighbor listをカットオフによって作成する
            cut_off: float
                カットオフ半径
            bond_length: list[list[float]]
                結合の長さ
        """
        count_bonds_lists = []
        for frame_idx in range(len(self.sf)):
            count_bonds_lists.append(
                self.sf[frame_idx].count_bonds(
                    mode=mode, cut_off=cut_off, bond_length=bond_length
                )
            )
        df_count_bonds = pd.DataFrame(count_bonds_lists).fillna(0).astype(int)
        df_count_bonds.index = self.get_step_nums()
        return df_count_bonds
