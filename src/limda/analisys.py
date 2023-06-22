from .SimulationFrames import SimulationFrames


#----------------------------------------------------------------------------------------------------
def split_sfs(sfs: SimulationFrames, list_size: int)->list[SimulationFrames]:
    """sfsを複数のsfsに分け, sfsのlistを返す。
        listのサイズを指定できる
    Parameters
    ----------
        sfs: SimulationFrames
            分割したいSimulationFrames
        list_size: int
            sfsを何分割するか
    Return val
    ----------
        sfs_list: list[SimulationFrames()]
        元のsfsを複数に分けたときのsfsから成るlist
    Example
    -------
        sfs  = SimulationFrames()  # len(sfs) = 10
        sfs_list = sfs.split_simulation_frames(3)
            ->sfs_listは3つのsfsからなるlistで、
                len(sfs_list[i]) = [4,3,3]
    """
    sfs_list = [SimulationFrames() for _ in range(list_size)]
    item_num_list = [int(len(sfs)/list_size) for _ in range(list_size)]
    for i in range(len(sfs) % list_size):
        item_num_list[i] += 1
    for idx, item in enumerate(item_num_list):
        for i in range(1, item+1):
            sfs_list[idx].sf.append(i+sum(item_num_list[0:idx]))
    for frames in sfs_list:
        frames.atom_symbol_to_type = sfs.atom_symbol_to_type
        frames.atom_type_to_symbol = sfs.atom_type_to_symbol
        frames.atom_type_to_mass = sfs.atom_type_to_mass
    return sfs_list
#-----------------------------------------------------------------------------------------------------------
def split_sfs_each(sfs: SimulationFrames, each_sfs_size: int, keep_remains:bool = False)->list[SimulationFrames]:
    """ sfsを複数のsfsに分け, sfsのlistを返す。
        sfs1つ1つサイズを指定できる。
    Parameters
    ----------
        sfs: SimulationFrames
            分割したいSimulationFrames
        each_sfs_size: int
            list内のsfsのサイズ
        keep_remains: bool
            残りを捨てるか、後ろにくっつけるか
    Return val
    ----------
        sfs_list: list[SimulationFrames()]
        元のsfsを複数に分けたときのsfsから成るlist
    Example
    -------
        sfs  = SimulationFrames()  # len(sfs) = 10
        sfs_list = sfs.split_simulation_frames(3)
            -> len(sfs_list[i]) = [3,3,3] (keep_remains = False)
               len(sfs_list[i]) = [3,3,3,1] (keep_remains = True)
    """
    list_size = int(len(sfs) / each_sfs_size)
    main_sfs = SimulationFrames()
    remain_sfs = SimulationFrames()
    remain = len(sfs)%each_sfs_size
    main_sfs.sf = sfs.sf[0:len(sfs)-remain]
    remain_sfs.sf = sfs.sf[len(sfs)-remain:len(sfs)]
    sfs_list = split_sfs(main_sfs, list_size)
    if keep_remains == True and remain != 0:
        remain_sfs.atom_symbol_to_type = sfs.atom_symbol_to_type
        remain_sfs.atom_type_to_symbol = sfs.atom_type_to_symbol
        remain_sfs.atom_type_to_mass = sfs.atom_type_to_mass
        sfs_list.append(remain_sfs)
    return sfs_list
