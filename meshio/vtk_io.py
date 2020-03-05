# -*- coding: utf-8 -*-
#
"""
I/O for VTK <https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf>.
"""
import logging
import numpy

from .mesh import Mesh
from .common import raw_from_cell_data


def __reshape_TENSOR_3D_FULL(value):
    v = value
    tens = numpy.array([[v[0], v[3], v[4]],
                        [v[4], v[1], v[5]],
                        [v[4], v[5], v[2]]])
    return tens


# https://www.vtk.org/doc/nightly/html/vtkCellType_8h_source.html
vtk_to_meshio_type = {
    0: "empty",
    1: "vertex",
    # 2: 'poly_vertex',
    3: "line",
    # 4: 'poly_line',
    5: "triangle",
    # 6: 'triangle_strip',
    # 7: 'polygon',
    # 8: 'pixel',
    9: "quad",
    10: "tetra",
    # 11: 'voxel',
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
    15: "penta_prism",
    16: "hexa_prism",
    21: "line3",
    22: "triangle6",
    23: "quad8",
    24: "tetra10",
    25: "hexahedron20",
    26: "wedge15",
    27: "pyramid13",
    28: "quad9",
    29: "hexahedron27",
    30: "quad6",
    31: "wedge12",
    32: "wedge18",
    33: "hexahedron24",
    34: "triangle7",
    35: "line4",
    #
    # 60: VTK_HIGHER_ORDER_EDGE,
    # 61: VTK_HIGHER_ORDER_TRIANGLE,
    # 62: VTK_HIGHER_ORDER_QUAD,
    # 63: VTK_HIGHER_ORDER_POLYGON,
    # 64: VTK_HIGHER_ORDER_TETRAHEDRON,
    # 65: VTK_HIGHER_ORDER_WEDGE,
    # 66: VTK_HIGHER_ORDER_PYRAMID,
    # 67: VTK_HIGHER_ORDER_HEXAHEDRON,
}
meshio_to_vtk_type = {v: k for k, v in vtk_to_meshio_type.items()}

vtk_type_to_numnodes = {
    0: 0,  # empty
    1: 1,  # vertex
    3: 2,  # line
    5: 3,  # triangle
    9: 4,  # quad
    10: 4,  # tetra
    12: 8,  # hexahedron
    13: 6,  # wedge
    14: 5,  # pyramid
    15: 10,  # penta_prism
    16: 12,  # hexa_prism
    21: 3,  # line3
    22: 6,  # triangle6
    23: 8,  # quad8
    24: 10,  # tetra10
    25: 20,  # hexahedron20
    26: 15,  # wedge15
    27: 13,  # pyramid13
    28: 9,  # quad9
    29: 27,  # hexahedron27
    30: 6,  # quad6
    31: 12,  # wedge12
    32: 18,  # wedge18
    33: 24,  # hexahedron24
    34: 7,  # triangle7
    35: 4,  # line4
}


# These are all VTK data types. One sometimes finds 'vtktypeint64', but
# this is ill-formed.
vtk_to_numpy_dtype_name = {
    "bit": "bool",
    "unsigned_char": "uint8",
    "char": "int8",
    "unsigned_short": "uint16",
    "short": "int16",
    "unsigned_int": "uint32",
    "int": "int32",
    "unsigned_long": "int64",
    "long": "int64",
    "float": "float32",
    "double": "float64",
}

numpy_to_vtk_dtype = {v: k for k, v in vtk_to_numpy_dtype_name.items()}


def read(filename):
    """Reads a VTK file."""
    with open(filename, "rb") as f:
        out = read_buffer(f)
    return out


def read_buffer(f):
    # initialize output data
    points = None
    field_data = {}
    cell_data_raw = {}
    point_data = {}

    # skip header and title
    f.readline()
    f.readline()

    data_type = f.readline().decode("utf-8").strip()
    assert data_type in ["ASCII", "BINARY"], "Unknown data type '{}'.".format(
        data_type
    )
    is_ascii = data_type == "ASCII"

    c = None
    ct = None

    # One of the problem in reading VTK files are POINT_DATA and CELL_DATA
    # fields. They can contain a number of SCALARS+LOOKUP_TABLE tables, without
    # giving and indication of how many there are. Hence, SCALARS must be
    # treated like a first-class section. To associate it with POINT/CELL_DATA,
    # we store the `active` section in this variable.
    active = None

    while True:
        line = f.readline().decode("utf-8")
        if not line:
            # EOF
            break

        line = line.strip()
        if len(line) == 0:
            continue

        split = line.split()
        section = split[0]

        if section == "DATASET":
            dataset_type = split[1]
            assert (
                dataset_type == "UNSTRUCTURED_GRID"
            ), "Only VTK UNSTRUCTURED_GRID supported (not {}).".format(
                dataset_type
            )

        elif section == "POINTS":
            active = "POINTS"
            num_points = int(split[1])
            data_type = split[2]
            points = _read_points(f, data_type, is_ascii, num_points)

        elif section == "CELLS":
            active = "CELLS"
            num_items = int(split[2])
            c = _read_cells(f, is_ascii, num_items)

        elif section == "CELL_TYPES":
            active = "CELL_TYPES"
            num_items = int(split[1])
            ct = _read_cell_types(f, is_ascii, num_items)

        elif section == "POINT_DATA":
            active = "POINT_DATA"
            num_items = int(split[1])

        elif section == "CELL_DATA":
            active = "CELL_DATA"
            num_items = int(split[1])

        elif section == "SCALARS":
            if active == "POINT_DATA":
                d = point_data
            else:
                assert active == "CELL_DATA", "Illegal SCALARS in section '{}'.".format(
                    active
                )
                d = cell_data_raw

            d.update(_read_scalar_field(f, num_items, split))

        elif section == "VECTORS":
            if active == "POINT_DATA":
                d = point_data
            else:
                assert active == "CELL_DATA", "Illegal SCALARS in section '{}'.".format(
                    active
                )
                d = cell_data_raw

            d.update(_read_vector_field(f, num_items, split, is_ascii))

        elif section == "TENSORS":
            if active == "POINT_DATA":
                d = point_data
            else:
                assert active == "CELL_DATA", "Illegal SCALARS in section '{}'.".format(
                    active
                )
                d = cell_data_raw

            d.update(_read_tensor_field(f, num_items, split, is_ascii))

        else:
            assert section == "FIELD", "Unknown section '{}'.".format(section)

            if active == "POINT_DATA":
                d = point_data
            else:
                assert active == "CELL_DATA", "Illegal FIELD in section '{}'.".format(
                    active
                )
                d = cell_data_raw

            d.update(_read_fields(f, int(split[2]), is_ascii))

    assert c is not None, "Required section CELLS not found."
    assert ct is not None, "Required section CELL_TYPES not found."

    cells, cell_data = translate_cells(c, ct, cell_data_raw)

    return Mesh(
        points, cells, point_data=point_data, cell_data=cell_data, field_data=field_data
    )


def _read_points(f, data_type, is_ascii, num_points):
    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    if is_ascii:
        #workaround for ABQ Error
        try:
            points = numpy.fromfile(f, count=num_points * 3, sep=" ", dtype=dtype)
        except:
            points = f.read().split()
            points = [float(i) for i in points]
            points = numpy.array(points)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        points = numpy.fromfile(f, count=num_points * 3, dtype=dtype)
        line = f.readline().decode("utf-8")
        assert line == "\n"

    return points.reshape((num_points, 3))


def _read_cells(f, is_ascii, num_items):
    if is_ascii:
        #workaround for ABQ Error
        try:
            c = numpy.fromfile(f, count=num_items, sep=" ", dtype=int)
        except:
            c = f.read().split()
            c = [int(i) for i in c]
            c = numpy.array(c)
            
    else:
        c = numpy.fromfile(f, count=num_items, dtype=">i4")
        line = f.readline().decode("utf-8")
        assert line == "\n"

    return c


def _read_cell_types(f, is_ascii, num_items):
    if is_ascii:
        #workaround for ABQ Error
        try:
            ct = numpy.fromfile(f, count=int(num_items), sep=" ", dtype=int)
        except:
            ct = f.read().split()
            ct = [int(i) for i in ct]
            ct = numpy.array(ct)
    else:
        # binary
        ct = numpy.fromfile(f, count=int(num_items), dtype=">i4")
        line = f.readline().decode("utf-8")
        # Sometimes, there's no newline at the end
        assert line.strip() == ""
    return ct


def _read_scalar_field(f, num_data, split):
    data_name = split[1]
    data_type = split[2]
    try:
        num_comp = int(split[3])
    except IndexError:
        num_comp = 1

    # The standard says:
    # > The parameter numComp must range between (1,4) inclusive; [...]
    assert 0 < num_comp < 5

    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    lt, _ = f.readline().decode("utf-8").split()
    assert lt == "LOOKUP_TABLE"
    #workaround for ABQ Error
    try:
        data = numpy.fromfile(f, count=num_data, sep=" ", dtype=dtype)
    except:
        data = f.read().split()
        data = [float(i) for i in data]
        data = numpy.array(data)
    return {data_name: data}


def _read_vector_field(f, num_data, split, is_ascii):
    data_name = split[1]
    data_type = split[2]

    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    if is_ascii:
        #workaround for ABQ Error
        try:
            data = numpy.fromfile(f, count=3 * num_data, sep=" ", dtype=dtype).reshape(-1, 3)
        except:
            data = f.read().split()
            data = [float(i) for i in data]
            data = numpy.array(data)
            data = data.reshape(-1,3)
    else:
        data = numpy.fromfile(f, count=num_data, sep="", dtype=dtype)

    return {data_name: data}


def _read_tensor_field(f, num_data, split, is_ascii):
    data_name = split[1]
    data_type = split[2]

    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    if is_ascii:
        #workaround for ABQ Error
        try:
            data = numpy.fromfile(f, count=9 * num_data, sep=" ", dtype=dtype).reshape(-1, 3, 3)
        except:
            data = f.read().split()
            data = [float(i) for i in data]
            data = numpy.array(data)
            data = data.reshape(-1,3,3)
    else:
        data = numpy.fromfile(f, count=9 * num_data, sep="", dtype=dtype).reshape(-1, 3, 3)

    return {data_name: data}


def _read_fields(f, num_fields, is_ascii):
    data = {}
    for _ in range(num_fields):
        name, shape0, shape1, data_type = f.readline().decode("utf-8").split()
        shape0 = int(shape0)
        shape1 = int(shape1)
        dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])

        if is_ascii:
            #workaround for ABQ Error
            try:
                dat = numpy.fromfile(f, count=shape0 * shape1, sep=" ", dtype=dtype)
            except:
                dat = f.read().split()
                dat = [float(i) for i in dat]
                dat = numpy.array(dat)
        else:
            # Binary data is big endian, see
            # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
            dtype = dtype.newbyteorder(">")
            dat = numpy.fromfile(f, count=shape0 * shape1, dtype=dtype)
            line = f.readline().decode("utf-8")
            assert line == "\n"

        if shape0 != 1:
            dat = dat.reshape((shape1, shape0))

        data[name] = dat

    return data


def translate_cells(data, types, cell_data_raw):
    # https://www.vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    # Translate it into the cells dictionary.
    # `data` is a one-dimensional vector with
    # (num_points0, p0, p1, ... ,pk, numpoints1, p10, p11, ..., p1k, ...

    # Collect types into bins.
    # See <https://stackoverflow.com/q/47310359/353337> for better
    # alternatives.
    bins = {u: numpy.where(types == u)[0] for u in numpy.unique(types)}

    # Deduct offsets from the cell types. This is much faster than manually
    # going through the data array. Slight disadvantage: This doesn't work for
    # cells with a custom number of points.
    numnodes = numpy.empty(len(types), dtype=int)
    for tpe, idx in bins.items():
        numnodes[idx] = vtk_type_to_numnodes[tpe]
    offsets = numpy.cumsum(numnodes + 1) - (numnodes + 1)
    assert numpy.all(numnodes == data[offsets])

    cells = {}
    cell_data = {}
    for tpe, b in bins.items():
        meshio_type = vtk_to_meshio_type[tpe]
        n = data[offsets[b[0]]]
        assert (data[offsets[b]] == n).all()
        indices = numpy.add.outer(offsets[b], numpy.arange(1, n + 1))
        cells[meshio_type] = data[indices]
        cell_data[meshio_type] = {key: value[b] for key, value in cell_data_raw.items()}

    return cells, cell_data


def write(filename, mesh, write_binary=True, encoding='utf-8'):
    if mesh.points.shape[1] == 2:
        logging.warning(
            "VTK requires 3D points, but 2D points given. "
            "Appending 0 third component."
        )
        mesh.points = numpy.column_stack(
            [mesh.points[:, 0], mesh.points[:, 1], numpy.zeros(mesh.points.shape[0])]
        )

    if not write_binary:
        logging.warning("VTK ASCII files are only meant for debugging.")

    if not write_binary:
        with open(filename, "w") as f:
            f.write("# vtk DataFile Version 4.2\n")
            f.write("{}\n".format(filename))
            #f.write("written by meshio v{}\n".format(__version__))
            f.write(("BINARY\n" if write_binary else "ASCII\n"))
            f.write("DATASET UNSTRUCTURED_GRID\n\n")

            # write points and cells
            _write_points(f, mesh.points, write_binary)
            _write_cells(f, mesh.cells, write_binary)

            # write point data
            if mesh.point_data:
                num_points = mesh.points.shape[0]
                f.write("POINT_DATA {}\n".format(num_points))
                _write_field_data(f, mesh.point_data, write_binary)

            # write cell data
            if mesh.cell_data:
                total_num_cells = sum([len(c) for c in mesh.cells.values()])
                cell_data_raw = raw_from_cell_data(mesh.cell_data)
                f.write("CELL_DATA {}\n".format(total_num_cells))
                _write_field_data(f, cell_data_raw, write_binary)

    elif write_binary:
        with open(filename, "wb") as f:
            f.write("# vtk DataFile Version 4.2\n".encode("utf-8"))
            f.write("{}\n".format(filename).encode("utf-8"))
            #f.write("written by meshio v{}\n".format(__version__))
            f.write(("BINARY\n" if write_binary else "ASCII\n").encode("utf-8"))
            f.write("DATASET UNSTRUCTURED_GRID\n\n".encode("utf-8"))

            # write points and cells
            _write_points(f, mesh.points, write_binary)
            _write_cells(f, mesh.cells, write_binary)

            # write point data
            if mesh.point_data:
                num_points = mesh.points.shape[0]
                f.write("POINT_DATA {}\n".format(num_points).encode("utf-8"))
                _write_field_data(f, mesh.point_data, write_binary)

            # write cell data
            if mesh.cell_data:
                total_num_cells = sum([len(c) for c in mesh.cells.values()])
                cell_data_raw = raw_from_cell_data(mesh.cell_data)
                f.write("CELL_DATA {}\n".format(total_num_cells).encode("utf-8"))
                _write_field_data(f, cell_data_raw, write_binary)


    return


def _write_points(f, points, write_binary):
    if write_binary:
        f.write(
            "POINTS {} {}\n".format(
                len(points), numpy_to_vtk_dtype[points.dtype.name]
            ).encode("utf-8")
        )
    else:
        f.write(
            "POINTS {} {}\n".format(
                len(points), numpy_to_vtk_dtype[points.dtype.name]
            )
        )

    if write_binary:
        # Binary data must be big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        points.astype(points.dtype.newbyteorder(">")).tofile(f, sep="")
    else:
        for point in points:
            #workaround to handle ABQ Error with tofile function
            try:
                point.tofile(f, sep=" ")
            except:
                f.write('%s %s %s ' %(point[0], point[1], point[2]))
            f.write('\n')
    if write_binary:
        f.write("\n".encode("utf-8"))
    else:
        f.write("\n")
    return


def _write_cells(f, cells, write_binary):
    total_num_cells = sum([len(c) for c in cells.values()])
    total_num_idx = sum([numpy.prod(c.shape) for c in cells.values()])
    # For each cell, the number of nodes is stored
    total_num_idx += total_num_cells
    if write_binary:
        f.write("CELLS {} {}\n".format(total_num_cells, total_num_idx).encode("utf-8"))
    else:
        f.write("CELLS {} {}\n".format(total_num_cells, total_num_idx))
    if write_binary:
        for key in cells:
            n = cells[key].shape[1]
            d = numpy.column_stack([numpy.ones(len(cells[key]))*n, cells[key]]).astype(
                numpy.dtype(">i4")
            )
            f.write(d.tostring())
        if write_binary:
            f.write("\n".encode("utf-8"))
    else:
        # ascii
        for key in cells:
            n = cells[key].shape[1]
            for cell in cells[key]:
                f.write(
                    (
                        " ".join(
                            ["{}".format(idx) for idx in numpy.concatenate([[n], cell])]
                        )
                        + "\n"
                    )
                )
    if write_binary:
        f.write("\n".encode("utf-8"))
    else:
        f.write("\n")
        f.write("\n")

    # write cell types
    if write_binary:
        f.write("CELL_TYPES {}\n".format(total_num_cells).encode("utf-8"))
    else:
        f.write("CELL_TYPES {}\n".format(total_num_cells))
    if write_binary:
        for key in cells:
            d = (numpy.ones(len(cells[key]))*meshio_to_vtk_type[key]).astype(
                numpy.dtype(">i4")
            )
            f.write(d.tostring())
        f.write("\n".encode("utf-8"))
    else:
        # ascii
        for key in cells:
            for _ in range(len(cells[key])):
                f.write("{}\n".format(meshio_to_vtk_type[key]))
        f.write("\n")
    return


def _write_field_data(f, data, write_binary):
    #f.write(("CELL DATA {}\n".format(len(data))))
    for name, values in data.items():
        if len(values.shape) == 1:
            num_tuples = values.shape[0]
            num_components = 1
            field_value_type = 'SCALARS'
        else:
            assert (
                len(values.shape) in [2,3]
            ), "Only one and two-dimensional field data supported."
            num_tuples = values.shape[0]
            if len(values.shape) == 2 and values.shape[-1] == 3:
                num_components = values.shape[1]
                field_value_type = 'VECTORS'
            elif len(values.shape) == 2 and values.shape[-1] == 6:
                field_value_type = 'TENSORS'
                values = numpy.array([__reshape_TENSOR_3D_FULL(x) for x in values])
                num_components = values.shape[1]*values.shape[2]
                num_dim = values.shape[1]
            else:
                num_components = values.shape[1]*values.shape[2]
                num_dim = values.shape[1]
                field_value_type = 'TENSORS'

        if " " in name:
            logging.warning(
                "VTK doesn't support spaces in field names. "
                "Renaming '%s' to '%s'.",
                name,
                name.replace(" ", "_"),
            )
            name = name.replace(" ", "_")
        if not write_binary:
            f.write(
                (
                    "{} {} {}\n".format(
                        field_value_type,
                        name,
                        numpy_to_vtk_dtype[values.dtype.name],
                    )
                )
            )
        else:
            f.write(
                (
                    "{} {} {}\n".format(
                        field_value_type,
                        name,
                        numpy_to_vtk_dtype[values.dtype.name],
                    )
                ).encode("utf-8")
            )
        if field_value_type == 'SCALARS':
                f.write(("LOOKUP_TABLE default\n").encode("utf-8"))
        if field_value_type == 'TENSORS':
            values = values.reshape(num_tuples*num_dim, num_dim)  # order
        if write_binary:
            values.astype(values.dtype.newbyteorder(">")).tofile(f, sep="")
        else:
            # ascii
            # values.tofile(f, sep=" ")
            for value in values:
                #workaround to handle ABQ Error with tofile function
                try:
                    value.tofile(f, sep=" ")
                except:
                    try:
                        for v in value:
                            f.write('%s ' %v)
                    except:
                        f.write('%s ' %value)
                f.write('\n')
            # numpy.savetxt(f, points)
    if write_binary:
        f.write("\n".encode("utf-8"))
    else:
        f.write("\n")
    return
