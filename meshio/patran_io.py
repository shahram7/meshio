# -*- coding: utf-8 -*-
#
"""
I/O for Patran files.
"""

import os
import numpy
from .mesh import Mesh

pat_to_meshio_type = {
        2: "line",
        3: "triangle",
        4: "quad",
        5: "tetra",
        7: "wedge",
        8: "hexa_prism",
        9: "pyramid",
}
meshio_to_pat_type = {v: k for k, v in pat_to_meshio_type.items()}


def read(filename):
    """Read a Patran *.pat file."""
    with open(filename, "r") as f:
        mesh, element_gids, point_gids = read_pat_buffer(f)

    # if *.ele file is present: Add cell data
    data_filename = filename.replace('.pat', '.ele')
    if os.path.isfile(data_filename):
        with open(data_filename, "r") as f:
            mesh = read_ele_buffer(f, mesh, element_gids)

    # if *.nod file is present: Add point data
    data_filename = filename.replace('.pat', '.nod')
    if os.path.isfile(data_filename):
        with open(data_filename, "r") as f:
            mesh = read_ele_buffer(f, mesh, point_gids)

    return mesh


def read_ele_buffer(f, mesh, element_gids):
    elem_id_map = {}
    for line, id in enumerate(element_gids):
        elem_id_map[id] = line

    name = f.readline().replace(' ', '_').rstrip('\n')
    dimensions = f.readline().split()
    N = int(dimensions[0])
    order = int(dimensions[-1])
    f.readline()

    array = numpy.zeros([N, order])

    for i in range(N):
        line = f.readline().split()
        ID = int(line[0])
        values = map(float, line[1:])
        line = elem_id_map[ID]
        array[line, :] = values

    shape_codes = mesh.cells.keys()
    arbitrary_shape_code = shape_codes[0]
    mesh.cell_data = {arbitrary_shape_code: {name: array}}
    return mesh


def read_nod_buffer(f, mesh, point_gids):
    node_id_map = {}
    for line, id in enumerate(point_gids):
        node_id_map[id] = line

    name = f.readline().replace(' ', '_').rstrip('\n')
    dimensions = f.readline().split()
    N = int(dimensions[0])
    order = int(dimensions[-1])
    f.readline()

    array = numpy.zeros([N, order])

    for i in range(N):
        line = f.readline().split()
        ID = int(line[0])
        values = map(float, line[1:])
        line = node_id_map[ID]
        array[line, :] = values

    mesh.point_data = {name: array}
    return mesh

def read_pat_buffer(f):
    # Initialize the optional data fields
    cells = {}
    points = []
    point_gids = []
    element_gids = []

    while True:
        line = f.readline()
        if not line:
            # EOF
            break

        card_data = line.split()
        if line.startswith(" 1"):
            points.append(_read_node(f))
            point_gids.append(int(card_data[1]))
        elif line.startswith(" 2"):
            lnodes = _read_cell(f)
            element_gids.append(int(card_data[1]))
            shape_code = int(card_data[2])
            key = pat_to_meshio_type[shape_code]
            if key in cells:
                cells[key] = numpy.vstack((cells[key], lnodes))
            else:
                cells[key] = lnodes

    points = numpy.array(points, dtype=float)
    point_gids = numpy.array(point_gids, dtype=int)

    cells = _scan_cells(point_gids, cells)

    return Mesh(points, cells), element_gids, point_gids


def _read_node(f):
    """ The node card contains the following:
    === ===== === ====== === ===
     1  ID    IV  KC
    === ===== === ====== === ===
    X   Y     Z
    ICF GTYPE NDF CONFIG CID PSP
    === ===== === ====== === ===
    """
    line = f.readline()
    entries = [line[i:i+16] for i in range(0, len(line), 16)]
    point = [float(coordinate) for coordinate in entries[:-1]]
    f.readline()
    return point


def _read_cell(f):
    """
    The element card contains the following:
    ====== ====== === ==== == == ==
    2      ID     IV  KC   N1 N2
    ====== ====== === ==== == == ==
    NODES  CONFIG PID CEID q1 q2 q3
    LNODES
    ADATA
    ====== ====== === ==== == == ==
    """
    f.readline()
    entries = f.readline().split()
    lnodes = map(int, entries)
    return numpy.array(lnodes)


def _scan_cells(point_gids, cells):
    for arr in cells.values():
        for value in numpy.nditer(arr, op_flags=["readwrite"]):
            value[...] = numpy.flatnonzero(point_gids == value)[0]
    return cells


def write(filename, mesh):
    # write a dummy
    with open(filename, "wt") as f:
        f.write("DUMMY")
