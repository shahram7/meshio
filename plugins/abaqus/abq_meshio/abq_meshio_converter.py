# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 18:04:22 2019

@author: wt4452
"""

from time import time

import numpy as np
import numpy.linalg as la

import meshio as mo

reload(mo)

in_abq = False
try:
    from abaqus import *
    from abaqusConstants import (
        NODAL,
        INTEGRATION_POINT,
        CENTROID,
        VECTOR,
        SCALAR,
        TENSOR_3D_FULL,
        TENSOR_3D_PLANAR,
        THREE_D,
        DEFORMABLE_BODY,
        TIME,
        MAGNITUDE,
        MISES,
        TRESCA,
        PRESS,
        INV3,
        MAX_PRINCIPAL,
    )
    from odbAccess import *
    from odbMaterial import *
    from odbSection import *

    in_abq = True
except ImportError:
    raise SystemError("Functions do only work in Abaqus")


def abaqus_to_meshio_type(element_type):
    """Map Abaqus elment type to meshio types.

    Parameters
    ----------
    element_type : str
        Abaqus element type (e.g C3D8R)

    Returns
    -------
    str
        Meshio element type (e.g. hexahedron)

    """

    # trusss
    if "T2D2" in element_type or "T3D2" in element_type:
        return "line"
    if "T2D3" in element_type or "T3D3" in element_type:
        return "line3"
    # beams
    if "B21" in element_type or "B31" in element_type:
        return "line"
    if "B22" in element_type or "B32" in element_type or "B33" in element_type:
        return "line3"
    # surfaces
    if "S4" in element_type or "R3D4" in element_type:
        return "quad"
    if "S8" in element_type:
        return "quad8"
    if "S8" in element_type:
        return "quad9"
    if "S3" in element_type or "M3D3" in element_type or "R3D3" in element_type:
        return "triangle"
    if "STRIA6" in element_type:
        return "triangle6"
    # volumes
    if "C3D8" in element_type or "EC3D8" in element_type or "SC8" in element_type:
        return "hexahedron"
    if "C3D20" in element_type:
        return "hexahedron20"
    if "C3D4" in element_type:
        return "tetra"
    if "C3D4H" in element_type:
        return "tetra4"
    if "C3D10" in element_type:
        return "tetra10"
    if "C3D6" in element_type:
        return "wedge"


meshio_to_abaqus_type = {
    "triangle": "S3R",
    "quad": "S4R",
    "hexahedron": "C3D8R",
    "tetra": "C3D4",
    "wedge": "C3D6",
}

# error messages
ERROR_NO_ODBObject = (
    "{} is no valid ODB object, please pass one of the "
    + "following: <odb>, <ODBAssembly>, <ODBInstance>"
)
ERROR_NO_MDBObject = (
    "{} is no valid MDB object, please pass one of the "
    + "following: <Part>, <PartInstance>, <Assembly>, "
    + "<Model>"
)
ERROR_NO_FIELD_DATA = "field output {} has no values in  " + "instance {}"
ERROR_ELSET_FIELD = (
    "field output {} is not defined on every element of " + "instance {}"
)
ERROR_RESHAPE_CELL_DATA = "{} cannot be reshaped into shape {}"
ERROR_DIFFERENT_ETYPES = (
    "different element types ({}) in {}. This feature " + "is not supported yet"
)


def __reshape_TENSOR_3D_FULL(value):
    v = value
    tens = np.array([[v[0], v[3], v[4]], [v[4], v[1], v[5]], [v[4], v[5], v[2]]])
    return tens


def __reshape_TENSOR_3D_PLANAR(value):
    v = value
    tens = np.array([[v[0], v[3], 0.0], [v[3], v[1], 0.0], [0.0, 0.0, v[2]]])
    return tens


# helper function for concatenating cell dictionaries
def __merge_numpy_dicts(dict1, dict2):
    new_dict = dict1.copy()
    old_keys = dict1.keys()
    new_keys = dict2.keys()

    for key in new_keys:
        if key in old_keys:
            new_dict[key] = np.vstack((dict1[key], dict2[key]))
        else:
            new_dict[key] = dict2[key]
    return new_dict


# helper function for concatenating cell_data dictionaries
def __merge_cellData_dicts(dict1, dict2):
    new_dict = dict1.copy()
    old_keys = dict1.keys()
    new_keys = dict2.keys()

    for key in new_keys:
        if key in old_keys:
            fO_1 = dict1[key]
            fO_2 = dict2[key]
            fO_1.update(fO_2)
        else:
            new_dict[key] = dict2[key]

    return new_dict


def __reshape_fieldOutputs(cell_data_field, allocation):
    """
    """
    new_cell_data_dict = {}
    new_field_name, field_names_to_reshape = allocation
    field_names_to_reshape = np.asarray(field_names_to_reshape)
    shape = field_names_to_reshape.shape

    field_names_to_reshape = field_names_to_reshape.flatten()

    assert shape in [(3, 3), (1, 3), (3, 1), (3,)], ERROR_RESHAPE_CELL_DATA.format(
        new_cell_data_dict.values(), shape
    )

    if not set(field_names_to_reshape).issubset(cell_data_field.keys()):
        return {}

    fields_to_reshape = np.array([cell_data_field[f] for f in field_names_to_reshape])
    fields_to_reshape = np.transpose(fields_to_reshape)
    n_values = fields_to_reshape.shape[0]
    if np.min(shape) > 1 and len(shape) > 1:
        new_field = fields_to_reshape.reshape((n_values, shape[0], shape[1]))
    else:
        new_field = fields_to_reshape

    new_cell_data_dict[new_field_name] = new_field

    return new_cell_data_dict


def convertMDBtoMeshio(mdbObject, **kwargs):
    """
    convertMDBtoMeshio(mdbObject, **kwargs)

    This function converts geometry information stored in Abaqus model database
    (mdb) to a meshio compatible representation.

    Parameters:
    ----------
    mdbObject : <'Part'> or <'PartInstance'> or <'Assembly'> or <'Model'>
        <'Part'> is defined on mdb.parts, whereas <'PartInstance'> is defined
        on the assembly level, <'Assembly'> may contain several part instance.
        When type <'Model'> is passed, its <'Assembly'> mebber is processed

    Returns
    -------
    Mesh : meshio Mesh object
        ready to write meshio Mesh objects
    """

    def convertInstance(mdbInstance, idx_shift=0):

        inst = mdbInstance
        nodes = inst.nodes
        elements = inst.elements

        n_nodes = len(nodes)

        # get node informations and coordinates
        node_labels, points = zip(*[(n.label, n.coordinates) for n in nodes])
        points = np.array(points)

        # create a lookup table to connect node labels and their array index
        nodeLU = {
            key: value
            for (key, value) in zip(node_labels, range(idx_shift, n_nodes + idx_shift))
        }

        # getting the elements is a bit more complex, since we have to sort by
        # type
        # firstly, we create an empty dict for storing the cell informations
        cells = {}
        cell_data = {}

        # loop over all elements
        for elem in elements:
            # get the connectivity
            con = [nodeLU[c + 1] for c in elem.connectivity]  # consider shift
            # get the type of element, convert to meshio representation
            etype = abaqus_to_meshio_type(str(elem.type))
            if etype in cells.keys():
                cells[etype].append(con)
                cell_data[etype]["ID"] = np.append(cell_data[etype]["ID"], elem.label)
            else:
                # create a new key for a new element set
                cells[etype] = [con]
                cell_data[etype] = {"ID": np.array([elem.label])}

        cells.update((key, np.array(cons)) for key, cons in cells.items())
        return points, cells, cell_data

    # if an Part or PartInstance is passed, call convertInstance once
    if str(type(mdbObject)) in ["<type 'Part'>", "<type 'PartInstance'>"]:
        points, cells, cell_data = convertInstance(mdbObject)

    # if an Assembly Object or a Model Object is passed, loop over
    # all instances
    elif str(type(mdbObject)) in ["<type 'Assembly'>", "<type 'Model'>"]:
        cells = {}
        cell_data = {}
        points = np.empty((0, 3))
        if str(type(mdbObject)) == "<type 'Model'>":
            rA = mdbObject.rootAssembly
        else:
            rA = mdbObject

        idx_shift = 0
        for inst_name in rA.instances.keys():
            inst = rA.instances[inst_name]
            points_, cells_, cell_data_ = convertInstance(inst, idx_shift)
            points = np.vstack((points, points_))
            cells = __merge_numpy_dicts(cells, cells_)
            cell_data = __merge_numpy_dicts(cell_data, cell_data_)
            idx_shift += len(points_)

    else:
        raise TypeError(ERROR_NO_MDBObject.format(mdbObject))

    return mo.Mesh(points, cells, cell_data=cell_data)


def convertODBtoMeshio(odbObject, frame, list_of_outputs=[], deformed=True, **kwargs):
    """
    convertODBtoMeshio(mdbObject, frame, list_of_outputs=None, **kwargs)

    This function converts geometry and result dat information stored in
    Abaqus Outbut database (odb) to a meshio compatible representation.

    Parameters:
    ----------
    odbObject : <'OdbInstance'> or <'OdbSet'> or <'OdbAssembly'> or <'odb'>
        <'OdbInstance'> is defined on odb.rootAssembly.instances, argument of
        type <'OdbSet'> must be an element Set as a member of an
        <'OdbInstance'> odject. <'OdbAssembly'> is the entire assembly
        containing several odbInstances and <'odb'> is the entire database
    frame: <'OdbFrame>'
        the frame containing the displacementField of the desire

    Returns
    -------
    Mesh : meshio Mesh object
        ready to write meshio Mesh objects
    """

    def convertInstance(odbInstance, frame, idx_shift=0, list_of_outputs=None):
        def processPointOutput(fO):

            # process node data
            print("processing " + fO.name)
            values = np.array([value.data for value in fO.values])
            if fO.type == SCALAR:
                point_data[fO.name] = values
            elif fO.type == VECTOR:
                point_data[fO.name] = values
            elif fO.type == TENSOR_3D_FULL:
                values_rs = np.array([__reshape_TENSOR_3D_FULL(v) for v in values])
                point_data[fO.name] = values_rs
            elif fO.type == TENSOR_3D_PLANAR:
                values_rs = np.array([__reshape_TENSOR_3D_PLANAR(v) for v in values])
                point_data[fO.name] = values_rs
            return

        def processCellOutput(fO):
            # process element data on several integration point
            if isElSet:
                fO_elem = fO.getSubset(region=eset)
            else:
                fO_elem = fO
            n_el_values = len(fO_elem.values)
            # check for availability of field output on each element
            assert n_el_values == n_elements, ERROR_ELSET_FIELD.format(
                field_name, inst_name
            )
            # use interpolation to output on centroid, to assert on result per
            # element
            fO_elem = fO_elem.getSubset(position=CENTROID)
            print("processing " + fO.name)
            etypes = set(
                [abaqus_to_meshio_type(etype) for etype in fO_elem.baseElementTypes]
            )
            # only one element type in fO
            assert len(etypes) == 1, ERROR_DIFFERENT_ETYPES.format(etypes, inst_name)
            if fO.type == TENSOR_3D_FULL:
                values_list = [
                    np.array(
                        [
                            __reshape_TENSOR_3D_FULL(value.data)
                            for value in fO_elem.values
                        ]
                    )
                ]
            elif fO.type == TENSOR_3D_PLANAR:
                values_list = [
                    np.array(
                        [
                            __reshape_TENSOR_3D_PLANAR(value.data)
                            for value in fO_elem.values
                        ]
                    )
                ]
            else:
                values_list = [np.array([value.data for value in fO_elem.values])]

            cell_data_ = {}
            for etype, values in zip(etypes, values_list):
                cell_data_[etype] = {fO.name: values}

            return cell_data_

        # assign
        inst = odbInstance
        inst_name = inst.name
        nodes = inst.nodes
        elements = inst.elements

        isElSet = str(type(odbInstance)) == "<type 'OdbSet'>" and hasattr(
            odbInstance, "elements"
        )

        if isElSet:
            eset = odbInstance
            elements = eset.elements
            inst_name = eset.name
            # get the instance that contains the element set
            # actually, we would have to check if *each* element in the set
            # is part of the same instance
            inst = eval(".".join(repr(elements[0]).split(".")[:-1]))
            nodes = inst.nodes

        # in the odb the initial coordinates are saved as a member to each
        # Node Object. i.o. to get the current coords, we have to extract the
        # displacement field and add the values
        dispField = frame.fieldOutputs["U"]
        dispField = dispField.getSubset(region=inst)
        disp_values = dispField.values
        nLU = inst.getNodeFromLabel  # built-in function

        n_nodes = len(nodes)
        n_elements = len(elements)
        # get node informations and coordinates
        node_labels, disp, x0 = zip(
            *[
                (value.nodeLabel, value.data, nLU(value.nodeLabel).coordinates)
                for value in disp_values
            ]
        )
        x0 = np.array(x0)
        disp = np.array(disp)
        if deformed:
            points = disp + x0
            print("Export deformed geometry.")
        else:
            points = x0
            print("Export undeformed geometry.")

        # create a lookup table to connect node labels and their array index
        nodeLU = {
            key: value
            for (key, value) in zip(node_labels, range(idx_shift, n_nodes + idx_shift))
        }

        # getting the elements is a bit more complex, since we have to sort by
        # type
        # firstly, we create an empty dict for storing the cell informations
        cells = {}

        # loop over all elements
        for elem in elements:
            # get the connectivity
            con = [nodeLU[c] for c in elem.connectivity]
            # get the type of element, convert to meshio representation
            etype = abaqus_to_meshio_type(str(elem.type))
            if etype in cells.keys():
                cells[etype].append(con)
            else:
                # create a new key for a new element set
                cells[etype] = [con]

        cells.update((key, np.array(cons)) for key, cons in cells.items())

        cell_data = {}
        point_data = {}

        # if field data is requested
        if list_of_outputs:
            for field_name in list_of_outputs:
                if type(field_name) == str and not field_name.lower().startswith(
                    "fdir"
                ):
                    fO = frame.fieldOutputs[field_name]
                    fO = fO.getSubset(region=inst)
                    n_values = len(fO.values)
                    if n_values > 0:
                        fO_location = fO.values[0].position  # NODAL, ELEMENT,
                        if fO_location == NODAL:
                            processPointOutput(fO)
                        elif fO_location in [CENTROID, INTEGRATION_POINT]:
                            cell_data_ = processCellOutput(fO)
                            cell_data = __merge_cellData_dicts(cell_data, cell_data_)
                    else:
                        print(ERROR_NO_FIELD_DATA.format(field_name, inst_name))

                elif type(field_name) in [list, set, tuple] and len(field_name) == 2:
                    new_field_names = field_name[1]
                    new_field_names = np.asarray(new_field_names)
                    new_field_names = new_field_names.flatten()
                    for fn in set(new_field_names):
                        fO = frame.fieldOutputs[fn]
                        fO = fO.getSubset(region=inst)
                        n_values = len(fO.values)
                        assert n_values > 0, ERROR_NO_FIELD_DATA.format(
                            field_name, inst_name
                        )

                        fO_location = fO.values[0].position  # NODAL, ELEMENT,
                        if fO_location == NODAL:
                            processPointOutput(fO)
                        elif fO_location in [CENTROID, INTEGRATION_POINT]:
                            cell_data_ = processCellOutput(fO)
                            cell_data = __merge_cellData_dicts(cell_data, cell_data_)

                    if fO_location == NODAL:
                        new_point_data = __reshape_fieldOutputs(point_data, field_name)
                        point_data.update(new_point_data)
                        for fn in set(new_field_names):
                            del point_data[fn]

                    if fO_location in [CENTROID, INTEGRATION_POINT]:
                        for etype, cd_dict in cell_data.items():
                            new_cd_dict = __reshape_fieldOutputs(cd_dict, field_name)
                            cd_dict.update(new_cd_dict)
                            for fn in set(new_field_names):
                                del cd_dict[fn]

        if "FDIR1" in list_of_outputs or "FDIR2" in list_of_outputs:
            # get initial fiber orientation from stress field
            stress = frame.fieldOutputs["S"].getSubset(region=eset)
            csys = np.asarray(stress.values[0].localCoordSystem)
            fdir1_0, fdir2_0 = csys[:2]

            def _computeDeformationGradient(con):
                """Compute the deformation gradient of the element."""
                assert (len(con)) in [3, 4], ""
                # coordinates in the initial configuration
                x0_coords = np.array([x0[c] for c in con])
                # coordinates in the current configuration
                x1_coords = np.array([points[c] for c in con])

                # compute the derivative of the iso-coordinates
                if len(con) == 3:  # linear triangle
                    B_xii = np.zeros((2, 3))
                    B_xii[0] = [-1.0, 1.0, 0.0]
                    B_xii[1] = [-1.0, 0.0, 1.0]
                else:  # linear quad at midpoint
                    B_xii = np.zeros((2, 4))
                    B_xii[0] = [-0.25, 0.25, 0.25, -0.25]
                    B_xii[1] = [-0.25, -0.25, 0.25, 0.25]
                # compute the Jacobians
                J_initial = np.dot(B_xii, x0_coords)
                J_initial_inv = np.dot(
                    la.inv(np.dot(J_initial, J_initial.T)), J_initial
                )
                J_current = np.dot(B_xii, x1_coords)
                # compute F as product of the Jacobians
                # F maps from inital to current configuration via
                # reference configuration.
                F = np.dot(J_current.T, J_initial_inv)
                return F

            for etype, cell_con in cells.items():

                def_grad = np.array(
                    [_computeDeformationGradient(con_idx) for con_idx in cell_con]
                )
                if "FDIR1" in list_of_outputs:
                    fdir1 = np.einsum("Ijk,k->Ij", def_grad, fdir1_0)
                    fdir1 = np.array([f_i / la.norm(f_i) for f_i in fdir1])
                    try:
                        cell_data[etype].update({"FDIR1": fdir1})
                    except KeyError:
                        cell_data[etype] = {"FDIR1": fdir1}
                if "FDIR2" in list_of_outputs:
                    fdir2 = np.einsum("Ijk,k->Ij", def_grad, fdir2_0)
                    fdir2 = np.array([f_i / la.norm(f_i) for f_i in fdir2])
                    try:
                        cell_data[etype].update({"FDIR2": fdir2})
                    except KeyError:
                        cell_data[etype] = {"FDIR2": fdir2}

                try:
                    cell_data[etype].update({"F": def_grad})
                except KeyError:
                    cell_data[etype] = {"F": def_grad}

        return points, cells, point_data, cell_data

    tic = time()
    if str(type(odbObject)) in ["<type 'OdbInstance'>", "<type 'OdbSet'>"]:
        odbInstance = odbObject
        points, cells, point_data, cell_data = convertInstance(
            odbInstance, frame, 0, list_of_outputs
        )

    elif str(odbObject.__class__) in ["<type 'Odb'>", "<type 'OdbAssembly'>"]:
        cells = {}
        points = np.empty((0, 3))

        point_data = {}
        cell_data = {}

        if str(type(odbObject)) == "<type 'Odb'>":
            rA = odbObject.rootAssembly
        else:
            rA = odbObject

        idx_shift = 0
        for inst_name in rA.instances.keys():
            inst = rA.instances[inst_name]

            if inst_name not in ["Assembly", "ASSEMBLY", "assembly"]:
                points_, cells_ = convertInstance(inst, frame, idx_shift)[:2]
                points = np.vstack((points, points_))
                cells = __merge_numpy_dicts(cells, cells_)
                idx_shift += len(points_)
    toc = time()
    print("took {} seconds".format(toc - tic))
    return mo.Mesh(points, cells, point_data, cell_data)


def convertMeshioToMDB(mesh, partname="test", modelname="Model-1", **kwargs):
    """
    This function creates a new part in the selected model database from
    the geometry information stored in a meshio Mesh object

    an OdbInstance object, defined in the 'rootAssembly' section
    (odb.rootAssembly.instances)
        * odbObject: an Abaqus OdbInstance object (<type 'OdbInstance'>)
    an elementSet of an OdbInstance object, defined in its ElementSet section
    (odbInstance.elementSets)
        * odbObject: an Abaqus OdbSet object (<type 'OdbSet'>)
    the entire assmebly of the model, containing several OdbInstances
        * odbObject: an Abaqus OdbAssembly Object (<type 'OdbAssembly'>)
    the entire Output Database
        * odbObject: an Abaqus OutputDataBase (<type 'Odb'>)

    Returns
    -------
    Mesh : meshio Mesh object
        ready to write meshio Mesh object
    """
    print(
        "Warning: only the geometry of the mesh can be converted to ABAQUS´"
        + "model database. Information on field data is lost."
    )
    assert type(mesh) == mo.mesh.Mesh, "No meshio Mesh instance"
    all_points = mesh.points
    all_cells = mesh.cells
    n_points = len(all_points)

    # convert nodes to abaqus comaptible representation
    nodeCoords = zip(all_points[:, 0], all_points[:, 1], all_points[:, 2])
    nodeLabels = list(range(1, n_points + 1))
    nodeData = [nodeLabels, nodeCoords]

    elementData = []
    element_label_shift = 1
    for etype, els in all_cells.items():
        # convert cells to abaqus compatible representation
        abq_element_type = meshio_to_abaqus_type[etype]
        element_labels = range(element_label_shift, len(els) + element_label_shift)

        # abaqus requires list and int data types
        # consider shifting node labels in connectivity
        element_con = [[int(i + 1) for i in x] for x in els]
        elementData.append([abq_element_type, element_labels, element_con])
        # continue counting in next iteration
        element_label_shift += len(els)

    model = mdb.models[modelname]
    partnames = model.parts.keys()
    i = 0
    while partname in partnames:
        if i == 0:
            print(
                "Warning: a part named {} is "
                / "already in model {}.".format(partname, modelname)
            )
        partname += "_{}".format(i)
        i += 1

    model.PartFromNodesAndElements(
        name=partname,
        dimensionality=THREE_D,
        type=DEFORMABLE_BODY,
        nodes=nodeData,
        elements=elementData,
    )


def convertMeshioToODB(mesh, odbname="test", filename="test.odb", **kwargs):

    all_points = mesh.points
    n_points = len(all_points)
    all_elements = mesh.cells
    all_node_data = mesh.point_data
    all_element_data = mesh.cell_data

    # creat a new odb
    odb = Odb(
        name=odbname,
        analysisTitle="ODB created from Meshio Instance",
        description="ODB created from Meshio Instance",
        path=filename,
    )

    # add section
    sCat = odb.SectionCategory(name="S5", description="Five-Layered Shell")

    # create part
    odb_part = odb.Part(name="part-1", embeddedSpace=THREE_D, type=DEFORMABLE_BODY)

    # get the nodes
    node_labels = range(1, 1 + n_points)
    node_list = zip(node_labels, *[all_points[:, i] for i in (0, 1, 2)])

    # add nodes to odb part
    odb_part.addNodes(nodeData=node_list, nodeSetName="nset-1")

    # get element data
    element_label_shift = 1
    element_labels_LU = {}
    for etype, els in all_elements.items():
        abq_element_type = meshio_to_abaqus_type[etype]
        element_labels = range(element_label_shift, len(els) + element_label_shift)
        # consider shifting node labels in connectivity
        element_con = [[int(i + 1) for i in x] for x in els]
        odb_part.addElements(
            labels=element_labels,
            connectivity=element_con,
            type=abq_element_type,
            elementSetName="{}".format(etype),
        )
        # continue counting in next iteration
        element_label_shift += len(els)
        element_labels_LU[etype] = element_labels

    # instance part
    odb_inst = odb.rootAssembly.Instance(name="part-1-1", object=odb_part)

    # create step and frame
    if all_node_data or all_element_data:
        dummy_step = odb.Step(
            name="step-1",
            domain=TIME,
            timePeriod=1.0,
            description="first analysis step",
        )
        dummy_frame = dummy_step.Frame(
            incrementNumber=0, frameValue=0.0, description="1st frame"
        )

    # write node data
    for nd_name, node_data in all_node_data.items():

        shape = node_data.shape
        if shape == (n_points,):
            field_type = SCALAR
        elif shape == (n_points, 3):
            field_type = VECTOR
        elif shape == (n_points, 3, 3):
            field_type = TENSOR_3D_FULL
        else:
            print("only scalar and vector data is supported atm")
            continue

        newField = dummy_frame.FieldOutput(
            name="{}".format(nd_name), description="{}".format(nd_name), type=field_type
        )
        if field_type == SCALAR:
            newField.addData(
                position=NODAL,
                instance=odb_inst,
                labels=node_labels,
                data=[[x] for x in node_data],
            )

        elif field_type == VECTOR:
            newField.setComponentLabels(("1", "2", "3"))
            newField.setValidInvariants((MAGNITUDE,))
            newField.addData(
                position=NODAL, instance=odb_inst, labels=node_labels, data=node_data
            )
        elif field_type == TENSOR_3D_FULL:
            # assume symmetry, keep only the relevant components,
            # reshape to abaqus order
            node_data = [x.flatten()[[0, 4, 8, 1, 2, 5]] for x in node_data]
            newField.setComponentLabels(("11", "22", "33", "12", "13", "23"))
            newField.setValidInvariants((MISES, TRESCA, PRESS, INV3, MAX_PRINCIPAL))
            newField.addData(
                position=NODAL,
                instance=odb_inst,
                labels=node_labels,
                data=[x.tolist() for x in node_data],
            )
        else:
            print("ERROR processing node output {}".format(nd_name))

    # write element data
    for etype, ed_dict in all_element_data.items():
        for ed_name, element_data in ed_dict.items():
            n_elements = len(element_data)
            shape = element_data.shape
            if shape == (n_elements,):
                field_type = SCALAR
            elif shape == (n_elements, 3):
                field_type = VECTOR
            elif shape == (n_elements, 3, 3):
                field_type = TENSOR_3D_FULL
            else:
                print("only scalar, vector and full3d data is supported atm")
                continue
            if ed_name not in dummy_frame.fieldOutputs.keys():
                # create new
                currentField = dummy_frame.FieldOutput(
                    name="{}".format(ed_name),
                    description="{}".format(ed_name),
                    type=field_type,
                )
            else:
                currentField = dummy_frame.fieldOutputs["{}".format(ed_name)]

            # add data to field_output
            if field_type == SCALAR:
                currentField.addData(
                    position=CENTROID,
                    instance=odb_inst,
                    labels=element_labels_LU[etype],
                    data=[[x] for x in element_data],
                )
            elif field_type == VECTOR:
                currentField.setComponentLabels(("1", "2", "3"))
                currentField.setValidInvariants((MAGNITUDE,))
                currentField.addData(
                    position=CENTROID,
                    instance=odb_inst,
                    labels=element_labels_LU[etype],
                    data=[x.tolist() for x in element_data],
                )
            elif field_type == TENSOR_3D_FULL:
                # assume symmetry, keep only the relevant components,
                # reshape to abaqus order
                element_data = [x.flatten()[[0, 4, 8, 1, 2, 5]] for x in element_data]
                currentField.setComponentLabels(("11", "22", "33", "12", "13", "23"))
                currentField.setValidInvariants(
                    (MISES, TRESCA, PRESS, INV3, MAX_PRINCIPAL)
                )
                currentField.addData(
                    position=CENTROID,
                    instance=odb_inst,
                    labels=element_labels_LU[etype],
                    data=[x.tolist() for x in element_data],
                )
            else:
                print("ERROR processing element output {}".format(ed_name))

    pathToODB = odb.path
    odb.save()
    odb.close()
    odb = openOdb(pathToODB)

    return odb
