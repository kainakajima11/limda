# distutils: language = c++

import numpy as np
import queue
from libc.stdlib cimport malloc
from libcpp.vector cimport vector
from libcpp.queue cimport queue
from libcpp cimport bool

cdef struct atom:
    # Cython上で原子のid, type, どこのmeshにいるか, 座標を記録する構造体
    int id
    int typ
    int mesh_id
    double pos[3]

cdef void make_catoms(vector[int] atoms_type, vector[vector[double]] atoms_pos, int atom_num, atom *catoms):
    # sfから持ってきた、原子のtype,positionをatom構造体として記録します。
    cdef int i
    for i in range(atom_num):
        catoms[i].typ = atoms_type[i]
        catoms[i].id = i
        catoms[i].pos[0] = atoms_pos[0][i]
        catoms[i].pos[1] = atoms_pos[1][i]
        catoms[i].pos[2] = atoms_pos[2][i]

cdef void make_mesh_size(vector[double] cell, double mesh_length, int mesh_size[3], double mesh_length_adjusted[3]):
    # x,y,zのmeshの個数を決めます。
    # またmeshの大きさが最適になるように調整します。
    cdef int i
    
    for i in range(3):
        mesh_size[i] = int(cell[i]/mesh_length)
        if mesh_size[i] < 3:
            mesh_size[i] = 3
        mesh_length_adjusted[i] = cell[i]/mesh_size[i]

cdef void make_mesh_id(int mesh_size[3], atom *catoms, double mesh_length_adjusted[3], int atom_num):
    # 原子がどこのmeshにいるのかを調べ、atom構造体のmesh_idに記録します。
    cdef:
        int mesh_num[3]
        int i,ax

    for i in range(atom_num):
        for ax in range(3):
            mesh_num[ax] = int(catoms[i].pos[ax] / mesh_length_adjusted[ax])
            if mesh_num[ax] == mesh_size[ax]:
                mesh_num[ax] -= 1

        catoms[i].mesh_id = mesh_num[2]*mesh_size[0]*mesh_size[1] + mesh_num[1]*mesh_size[0] + mesh_num[0]

cdef vector[vector[int]] get_append_mesh(atom *catoms, vector[vector[int]] append_mesh, int atom_num):
    # meshそれぞれにどの原子がいるかをappend_meshに記録します.
    cdef int i
    for i in range(atom_num):
        append_mesh[catoms[i].mesh_id].push_back(catoms[i].id)

    return append_mesh

cdef vector[vector[int]] search_neighbors(atom *catoms,
                                          vector[vector[int]] append_mesh,
                                          int mesh_size[3],
                                          vector[vector[int]] neighbor_list,
                                          vector[vector[double]] bond_length,
                                          vector[double] cell):
    # 近接meshを探索して、原子の結合listを作成します。
    # bond_lengthを使用します.
    cdef:
        double dx[3]
        int own_mesh_len, serche_mesh_len
        atom own, search
        int own_mesh_id, search_mesh_id, iid, jid, dim

    add_idxes = (
        [0,0,0],
        [1,0,0],
        [-1,1,0],
        [0,1,0],
        [1,1,0],
        [-1,0,1],
        [0,0,1],
        [1,0,1],
        [-1,-1,1],
        [0,-1,1],
        [1,-1,1],
        [-1,1,1],
        [0,1,1],
        [1,1,1],
    )
    for own_mesh_id in range(mesh_size[0]*mesh_size[1]*mesh_size[2]):
        own_mesh_len = len(append_mesh[own_mesh_id])
        for add_idx in add_idxes:
            search_mesh_id = (own_mesh_id%mesh_size[0] + add_idx[0])%mesh_size[0] \
                + mesh_size[0]*((own_mesh_id/mesh_size[0] + add_idx[1])%mesh_size[1]) \
                + mesh_size[0]*mesh_size[1]*((own_mesh_id/(mesh_size[0]*mesh_size[1]) + add_idx[2])%mesh_size[2])
            search_mesh_len = len(append_mesh[search_mesh_id])
            for iid in range(own_mesh_len):
                own = catoms[append_mesh[own_mesh_id][iid]]
                start = (search_mesh_id == own_mesh_id)
                for jid in range(start*(iid+1), search_mesh_len):
                    search = catoms[append_mesh[search_mesh_id][jid]]
                    for dim in range(3):
                        dx[dim] = search.pos[dim] - own.pos[dim]
                        if dx[dim] < -cell[dim]/2:
                            dx[dim] += cell[dim]
                        if cell[dim]/2 < dx[dim]:
                            dx[dim] -= cell[dim]
                    if dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2] <= bond_length[own.typ-1][search.typ-1]*bond_length[own.typ-1][search.typ-1]:
                        neighbor_list[own.id].push_back(search.id)
                        neighbor_list[search.id].push_back(own.id)
    return neighbor_list

cdef vector[vector[int]] make_neighbor_list(vector[int] atoms_type,
                                            vector[vector[double]] atoms_pos,
                                            double mesh_length,
                                            int atom_num,
                                            vector[vector[double]] bond_length,
                                            vector[double] cell):
    cdef:
        atom *catoms = <atom *> malloc(atom_num * sizeof(atom))
        int mesh_size[3]
        double mesh_length_adjusted[3]
        vector[vector[int]] append_mesh
        vector[vector[int]] neighbor_list

    make_catoms(atoms_type, atoms_pos, atom_num, catoms)
    make_mesh_size(cell, mesh_length, mesh_size, mesh_length_adjusted)
    make_mesh_id(mesh_size, catoms, mesh_length_adjusted, atom_num)
    append_mesh.resize(mesh_size[0]*mesh_size[1]*mesh_size[2])
    append_mesh = get_append_mesh(catoms, append_mesh, atom_num)
    neighbor_list.resize(atom_num)
    neighbor_list = search_neighbors(catoms, append_mesh, mesh_size, neighbor_list, bond_length, cell)
    return neighbor_list

def get_neighbor_list_using_cython(vector[int] atoms_type,
                                   vector[vector[double]] atoms_pos,
                                   double mesh_length,
                                   int atom_num,
                                   vector[vector[double]] bond_length,
                                   vector[double] cell):
    return make_neighbor_list(atoms_type,atoms_pos,mesh_length,atom_num,bond_length,cell)
