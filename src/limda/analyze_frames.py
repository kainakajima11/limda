import pandas as pd


class AnalyzeFrames:
    def __init__(self):
        pass

    def count_mols(
        self,
        mode: str = "bond_length",
        cut_off: float = None,
        bond_length: list[list[float]] = None,
    ) -> pd.DataFrame:
        """分子数を数える
以下は0 ~ 8000000 stepの分子を数えた例
```
         C4H11O4  C2H5O2  C2H7O2  C2H4O2  H2O1  C4H10O4  C2H6O1  C2H4  C2H5O1
0              4      65       2       9     1        1       1     1       1  
1000000        0      28       0      11     3        0      10     3       6
2000000        0      16       0      10     4        1      15     6       4
3000000        0      16       0       8     7        1      18     9       9
4000000        0      19       1      14     6        0      19    14      14
5000000        1      10       0      13     6        1      25    19      12
6000000        0       5       0      13     5        1      25    18      14
7000000        0       8       0      10     7        1      31    18      11
8000000        0       5       0       7     6        1      40    18       7
```
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
        columns = list(df_count_mols.columns)
        columns.sort(key= lambda col: df_count_mols[col].max(), reverse=True) # 出現率の多い分子から表示
        df_count_mols = df_count_mols[columns]
        return df_count_mols

    def count_bonds(
        self,
        mode: str = "bond_length",
        cut_off: float = None,
        bond_length: list[list[float]] = None,
    ) -> pd.DataFrame:
        """結合数を数える
以下は0 ~ 8000000 stepの結合数を数えた例
```
          C-C   C-H   C-O  C-N  C-Si  H-H  H-O  H-N  H-Si  O-O  O-N  O-Si
0        1000  3961  1988    0     0    0  327    0     0    0    0     0  
1000000  1008  3931  1915    0     0    0  603    0     0    0    0     0
2000000  1015  3904  1882    0     0    0  611    0     0    0    0     0
3000000  1022  3875  1824    0     0    0  619    0     0    0    0     0
4000000  1036  3825  1757    0     0    0  575    0     0    0    0     0
5000000  1047  3779  1691    0     0    0  577    0     0    0    0     0
6000000  1054  3745  1660    0     0    0  551    0     0    0    0     0
7000000  1064  3718  1617    0     0    0  540    0     0    0    0     0
8000000  1076  3675  1572    0     0    0  514    0     0    0    0     0

         N-N   N-Si  Si-Si
0        129  29063    915
1000000  123  28943    888
2000000  126  28875    874
3000000  120  28802    859
4000000  127  28636    847
5000000  121  28513    847
6000000  126  28435    856
7000000  118  28342    848
8000000  120  28255    818
```
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
