# -*- coding: utf-8 -*-
#
"""I/O for Patran files.

The geometry is read from the provided *.pat file. If a file with the same
name and ending *.ele or *.nod is provided, cell data and point data is added.

This feature is designed to convert Patran files generated from Modlflow(See
https://knowledge.autodesk.com/support/moldflow-insight/learn-explore/
caas/CloudHelp/cloudhelp/2018/ENU/MoldflowInsight/files/
GUID-BCC20E1A-12EA-428F-95F5-C1E4BC1E416C-htm.html and
https://knowledge.autodesk.com/support/moldflow-insight/learn-explore/
caas/sfdcarticles/sfdcarticles/How-to-export-fiber-orientation-results
-in-XML-or-Patran-format-from-Moldflow.html)

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


def read(filename, ele_filename=None, nod_filename=None, scale=1.0):
    """Read a Patran *.pat file.

    If a *.ele file or *.nod file is provided or if it has the same name as
    the *.pat file, these files are used to fill data fields. Such files are
    exported by Modlflow for example.

    Args
    ----
        filename (str): Patran filename that should be read

        ele_filename (str, optional): element-wise data file

        nod_filename (str, optional): node-wise data file.

        scale (float): scale factor for nodes
    """
    with open(filename, "r") as f:
        mesh, element_gids, point_gids = read_pat_buffer(f, scale)

    # if *.ele file is present: Add cell data
    ele_filename = ele_filename or filename.replace('.pat', '.ele')
    if os.path.isfile(ele_filename):
        with open(ele_filename, "r") as f:
            mesh = read_ele_buffer(f, mesh, element_gids)

    # if *.nod file is present: Add point data
    nod_filename = nod_filename or filename.replace('.pat', '.nod')
    if os.path.isfile(nod_filename):
        with open(nod_filename, "r") as f:
            mesh = read_nod_buffer(f, mesh, point_gids)

    return mesh


def read_ele_buffer(f, mesh, element_gids):
    """Read element based data file."""

    name = f.readline().replace(' ', '_').rstrip('\n')
    dimensions = f.readline().split()
    N = int(dimensions[0])
    f.readline()

    data = {}

    for i in range(N):
        line = f.readline().split()
        ID = int(line[0])
        values = map(float, line[1:])
        data[ID] = numpy.array(values)

    for elem_type in mesh.cells.keys():
        for gid in element_gids[elem_type]:
            try:
                values = data[gid]
                if elem_type in mesh.cell_data.keys():
                    mesh.cell_data[elem_type][name].append(values)
                else:
                    mesh.cell_data[elem_type] = {name: [values]}
            except KeyError:
                print("No data attached to Element %d" % gid)

    return mesh


def read_nod_buffer(f, mesh, point_gids):
    """Read node based data file."""
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


def read_pat_buffer(f, scale):
    # Initialize the optional data fields
    cells = {}
    points = []
    point_gids = []
    element_gids = {}

    while True:
        line = f.readline()
        if not line:
            # EOF
            break

        card_data = line.split()
        if line.startswith(" 1"):
            points.append(_read_node(f, scale))
            point_gids.append(int(card_data[1]))
        elif line.startswith(" 2"):
            lnodes = _read_cell(f)
            shape_code = int(card_data[2])
            key = pat_to_meshio_type[shape_code]
            if key in cells:
                cells[key] = numpy.vstack((cells[key], lnodes))
                element_gids[key].append(int(card_data[1]))
            else:
                cells[key] = lnodes
                element_gids[key] = [int(card_data[1])]
        elif line.startswith(" 4"):
            # do not read cross section properties.
            junk = f.readline()

    points = numpy.array(points, dtype=float)
    point_gids = numpy.array(point_gids, dtype=int)

    cells = _scan_cells(point_gids, cells)

    return Mesh(points, cells), element_gids, point_gids


def _read_node(f, scale):
    """Read a node card.

    The node card contains the following:
    === ===== === ====== === ===
     1  ID    IV  KC
    === ===== === ====== === ===
    X   Y     Z
    ICF GTYPE NDF CONFIG CID PSP
    === ===== === ====== === ===
    """
    line = f.readline()
    entries = [line[i:i+16] for i in range(0, len(line), 16)]
    point = [scale*float(coordinate) for coordinate in entries[:-1]]
    f.readline()
    return point


def _read_cell(f):
    """Read a cell card.

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
    lnodes = numpy.trim_zeros(lnodes)
    return numpy.array(lnodes)


def _scan_cells(point_gids, cells):
    for arr in cells.values():
        for value in numpy.nditer(arr, op_flags=["readwrite"]):
            value[...] = numpy.flatnonzero(point_gids == value)[0]
    return cells


def write(filename, mesh):
    """Write a dummy for now."""
    with open(filename, "wt") as f:
        f.write("DUMMY")
