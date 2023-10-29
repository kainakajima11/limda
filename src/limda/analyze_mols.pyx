# distutils: language = c++

import numpy as np
import queue
from libc.stdlib cimport malloc
from libcpp.vector cimport vector
from libcpp.queue cimport queue
from libcpp cimport bool

cdef vector[vector[int]] get_mols_list(vector[vector[int]] neighbor_list, int atom_num):
    cdef:
        vector[bool] visited
        queue[int] que
        int start_atom_idx
        vector[vector[int]] mols_list
        int mol_num

    mol_num = 0
    visited.resize(atom_num)

    for start_atom_idx in range(atom_num):
        if visited[start_atom_idx]:
            continue
        mol_num += 1
        mols_list.push_back([])
        que.push(start_atom_idx)
        while not que.empty():
            now = que.front()
            que.pop()
            if visited[now]:
                continue
            visited[now] = True
            mols_list[mol_num - 1].push_back(now)
            for nex in neighbor_list[now]:
                if visited[nex]:
                    continue
                que.push(nex)
    return mols_list
#------------------------------------------------------------------------------------------------------------------
def get_mols_list_using_cython(vector[vector[int]] neighbor_list, int atom_num):
    return get_mols_list(neighbor_list, atom_num)
#------------------------------------------------------------------------------------------------------------------