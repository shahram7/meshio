# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 18:04:22 2019

@author: wt4452
"""

import numpy as np
import meshio as mo

reload(mo)

in_abq = False
try:
    from abaqus import *
    from abaqusConstants import (NODAL, INTEGRATION_POINT, CENTROID,
                                 VECTOR, SCALAR, TENSOR_3D_FULL,
                                 TENSOR_3D_SURFACE, TENSOR_3D_PLANAR,
                                 THREE_D, DEFORMABLE_BODY, LONG_TERM,
                                 ISOTROPIC, TIME, MAGNITUDE, MISES, TRESCA,
                                 PRESS, INV3, MAX_PRINCIPAL)
    from odbAccess import *
    from odbMaterial import *
    from odbSection import *
    in_abq = True
except ImportError:
    raise SystemError('Functions do only work in Abaqus')

# not complete, has to be completed with new elements (membrane, rigid, ...)
abaqus_to_meshio_type = {
    # trusss
    "T2D2": "line",
    "T2D2H": "line",
    "T2D3": "line3",
    "T2D3H": "line3",
    "T3D2": "line",
    "T3D2H": "line",
    "T3D3": "line3",
    "T3D3H": "line3",
    # beams
    "B21": "line",
    "B21H": "line",
    "B22": "line3",
    "B22H": "line3",
    "B31": "line",
    "B31H": "line",
    "B32": "line3",
    "B32H": "line3",
    "B33": "line3",
    "B33H": "line3",
    # surfaces
    "S4": "quad",
    "S4R": "quad",
    "S4RS": "quad",
    "S4RSW": "quad",
    "S4R5": "quad",
    "R3D4": "quad",
    "S8R": "quad8",
    "S8R5": "quad8",
    "S9R5": "quad9",
    # "QUAD": "quad",
    # "QUAD4": "quad",
    # "QUAD5": "quad5",
    # "QUAD8": "quad8",
    # "QUAD9": "quad9",
    #
    "STRI3": "triangle",
    "S3": "triangle",
    "S3R": "triangle",
    "S3RS": "triangle",
    "M3D3": "triangle",
    "R3D3": "triangle",
    # "TRI7": "triangle7",
    # 'TRISHELL': 'triangle',
    # 'TRISHELL3': 'triangle',
    # 'TRISHELL7': 'triangle',
    #
    "STRI65": "triangle6",
    # 'TRISHELL6': 'triangle6',
    # volumes
    "C3D8": "hexahedron",
    "C3D8H": "hexahedron",
    "C3D8I": "hexahedron",
    "C3D8IH": "hexahedron",
    "C3D8R": "hexahedron",
    "C3D8RH": "hexahedron",
    # "HEX9": "hexahedron9",
    "C3D20": "hexahedron20",
    "C3D20H": "hexahedron20",
    "C3D20R": "hexahedron20",
    "C3D20RH": "hexahedron20",
    # "HEX27": "hexahedron27",
    #
    "C3D4": "tetra",
    "C3D4H": "tetra4",
    # "TETRA8": "tetra8",
    "C3D10": "tetra10",
    "C3D10H": "tetra10",
    "C3D10I": "tetra10",
    "C3D10M": "tetra10",
    "C3D10MH": "tetra10",
    # "TETRA14": "tetra14",
    #
    # "PYRAMID": "pyramid",
    "C3D6": "wedge",
}

# not complete --> update soon
meshio_to_abaqus_type = {
    'triangle': 'S3R',
    'quad': 'S4R',
    'hexahedron': 'C3D8R',
    'tetra': 'C3D4',
    'wedge': 'C3D6',

}

# error messages
ERROR_NO_ODBObject = '{} is no valid ODB object, please pass one of the ' \
                     + 'following: <odb>, <ODBAssembly>, <ODBInstance>'
ERROR_NO_MDBObject = '{} is no valid MDB object, please pass one of the ' \
                     + 'following: <Part>, <PartInstance>, <Assembly>, ' \
                     + '<Model>'
ERROR_NO_FIELD_DATA = 'field output {} has no values in  ' \
                    + 'instance {}'
ERROR_ELSET_FIELD = 'field output {} is not defined on every element of ' \
                    + 'instance {}'
ERROR_RESHAPE_CELL_DATA = '{} cannot be reshaped into shape {}'
ERROR_DIFFERENT_ETYPES = 'different element types ({}) in {}. This feature ' \
                         + 'is not supported yet'


def __reshape_TENSOR_3D_FULL(value):
    v = value
    tens = np.array([[v[0], v[3], v[4]],
                     [v[4], v[1], v[5]],
                     [v[4], v[5], v[2]]])
    return tens


def __reshape_TENSOR_3D_PLANAR(value):
    v = value
    tens = np.array([[v[0], v[3],   0.],
                     [v[3], v[1],   0.],
                     [0.,   0.,   v[2]]])
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
    '''
    '''
    new_cell_data_dict = {}
    new_field_name, field_names_to_reshape = allocation
    field_names_to_reshape = np.asarray(field_names_to_reshape)
    shape = field_names_to_reshape.shape

    field_names_to_reshape = field_names_to_reshape.flatten()

    assert shape in [(3, 3), (1, 3), (3, 1), (3,)], \
        ERROR_RESHAPE_CELL_DATA.format(new_cell_data_dict.values(), shape)

    if not set(field_names_to_reshape).issubset(cell_data_field.keys()):
        return {}

    fields_to_reshape = np.array([cell_data_field[f]
                                  for f in field_names_to_reshape])
    fields_to_reshape = np.transpose(fields_to_reshape)
    n_values = fields_to_reshape.shape[0]
    if np.min(shape) > 1 and len(shape) > 1:
        new_field = fields_to_reshape.reshape((n_values,
                                               shape[0], shape[1]))
    else:
        new_field = fields_to_reshape

    new_cell_data_dict[new_field_name] = new_field

    return new_cell_data_dict


def convertMDBtoMeshio(mdbObject, **kwargs):
    '''
    This function converts geometry information stored in Abaqus model database
    (mdb) to a meshio compatible representation

    Three different inputs can be passed and processed:

    a part object, defined in the 'Part' module (mdb.parts)
        * mdbObject: an Abaqus Part object (<type 'Part'>)
    an instance of an part object, defined in the 'Assembly' module
    (mdb.rootAssembly.instances)
        * mdbObject: an Abaqus PartInstance object (<type 'PartInstance'>)
    the entire assmebly of the model, containing several PartInstances
        * mdbObject: an Abaqus Assembly Object (<type 'Assembly'>)

    Returns
    -------
    Mesh : meshio Mesh object
        ready to write meshio Mesh object
    '''

    def convertInstance(mdbInstance, idx_shift=0):

        inst = mdbInstance
        nodes = inst.nodes
        elements = inst.elements

        n_nodes = len(nodes)

        # get node informations and coordinates
        node_labels, points = zip(*[(n.label, n.coordinates) for n in nodes])
        points = np.array(points)

        # create a lookup table to connect node labels and their array index
        nodeLU = {key: value for (key, value) in zip(node_labels,
                  range(idx_shift, n_nodes+idx_shift))}

        # getting the elements is a bit more complex, since we have to sort by
        # type
        # firstly, we create an empty dict for storing the cell informations
        cells = {}

        # loop over all elements
        for elem in elements:
            # get the connectivity
            con = [nodeLU[c+1] for c in elem.connectivity]  # consider shift
            # get the type of element, convert to meshio representation
            etype = abaqus_to_meshio_type[str(elem.type)]
            if etype in cells.keys():
                cells[etype].append(con)
            else:
                # create a new key for a new element set
                cells[etype] = [con]

        cells.update((key, np.array(cons)) for key, cons in cells.items())
        return points, cells

    # if an Part or PartInstance is passed, call convertInstance once
    if str(type(mdbObject)) in ["<type 'Part'>", "<type 'PartInstance'>"]:
        points, cells = convertInstance(mdbObject)

    # if an Assembly Object or a Model Object is passed, loop over
    # all instances
    elif str(type(mdbObject)) in ["<type 'Assembly'>", "<type 'Model'>"]:
        cells = {}
        points = np.empty((0, 3))
        if str(type(mdbObject)) == "<type 'Model'>":
            rA = mdbObject.rootAssembly
        else:
            rA = mdbObject

        idx_shift = 0
        for inst_name in rA.instances.keys():
            inst = rA.instances[inst_name]
            points_, cells_ = convertInstance(inst, idx_shift)
            points = np.vstack((points, points_))
            cells = __merge_numpy_dicts(cells, cells_)
            idx_shift += len(points_)

    else:
        raise TypeError(ERROR_NO_MDBObject.format(mdbObject))

    return mo.Mesh(points, cells)


def convertODBtoMeshio(odbObject, frame, list_of_outputs=None, **kwargs):

    def convertInstance(odbInstance, frame, idx_shift=0,
                        list_of_outputs=None):

        def processPointOutput(fO):

            # process node data
            print 'processing ' + fO.name
            values = np.array([value.data for value in fO.values])
            if fO.type == SCALAR:
                point_data[fO.name] = values
            elif fO.type == VECTOR:
                point_data[fO.name] = values
            elif fO.type == TENSOR_3D_FULL:
                values_rs = np.array([__reshape_TENSOR_3D_FULL(v)
                                      for v in values])
                point_data[fO.name] = values_rs
            elif fO.type == TENSOR_3D_PLANAR:
                values_rs = np.array([__reshape_TENSOR_3D_PLANAR(v)
                                      for v in values])
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
            assert n_el_values == n_elements, (ERROR_ELSET_FIELD
                                               .format(field_name,
                                                       inst_name))
            # use interpolation to output on centroid, to assert on result per
            # element
            fO_elem = fO_elem.getSubset(position=CENTROID)
            print 'processing ' + fO.name
            etypes = set([abaqus_to_meshio_type[etype]
                          for etype in fO_elem.baseElementTypes])
            # only one element type in fO
            assert len(etypes) == 1,  ERROR_DIFFERENT_ETYPES.format(etypes,
                                                                    inst_name)
            if fO.type == TENSOR_3D_FULL:
                values_list = [np.array([__reshape_TENSOR_3D_FULL(value.data)
                                         for value in fO_elem.values])]
            elif fO.type == TENSOR_3D_PLANAR:
                values_list = [np.array([__reshape_TENSOR_3D_PLANAR(value.data)
                               for value in fO_elem.values])]
            else:
                values_list = [np.array([value.data
                               for value in fO_elem.values])]

            cell_data_ = {}
            for etype, values in zip(etypes, values_list):
                cell_data_[etype] = {fO.name: values}

            return cell_data_

        # assign
        inst = odbInstance
        inst_name = inst.name
        nodes = inst.nodes
        elements = inst.elements

        isElSet = str(type(odbInstance)) == "<type 'OdbSet'>" and \
                                            hasattr(odbInstance, 'elements')

        if isElSet:
            eset = odbInstance
            elements = eset.elements
            inst_name = eset.name
            # get the instance that contains the element set
            # actually, we would have to check if *each* element in the set
            # is part of the same instance
            inst = eval('.'.join(repr(elements[0]).split('.')[:-1]))
            nodes = inst.nodes

        # in the odb the initial coordinates are saved as a member to each
        # Node Object. i.o. to get the current coords, we have to extract the
        # displacement field and add the values
        dispField = frame.fieldOutputs['U']
        dispField = dispField.getSubset(region=inst)
        disp_values = dispField.values
        nLU = inst.getNodeFromLabel  # built-in function

        n_nodes = len(nodes)
        n_elements = len(elements)
        # get node informations and coordinates
        node_labels, points = zip(*[(value.nodeLabel, value.data +
                                     nLU(value.nodeLabel).coordinates)
                                  for value in disp_values])

        points = np.array(points)

        # create a lookup table to connect node labels and their array index
        nodeLU = {key: value for (key, value) in zip(node_labels,
                  range(idx_shift, n_nodes+idx_shift))}

        # getting the elements is a bit more complex, since we have to sort by
        # type
        # firstly, we create an empty dict for storing the cell informations
        cells = {}

        # loop over all elements
        for elem in elements:
            # get the connectivity
            con = [nodeLU[c] for c in elem.connectivity]
            # get the type of element, convert to meshio representation
            etype = abaqus_to_meshio_type[str(elem.type)]
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
                if type(field_name) == str:
                    fO = frame.fieldOutputs[field_name]
                    fO = fO.getSubset(region=inst)
                    n_values = len(fO.values)
                    assert n_values > 0, ERROR_NO_FIELD_DATA.format(field_name,
                                                                    inst_name)
                    fO_location = fO.values[0].position  # NODAL, ELEMENT, ...
                    if fO_location == NODAL:
                        processPointOutput(fO)
                    elif fO_location in [CENTROID, INTEGRATION_POINT]:
                        cell_data_ = processCellOutput(fO)
                        cell_data = __merge_cellData_dicts(cell_data,
                                                           cell_data_)

                elif type(field_name) in [list, set,
                                          tuple] and len(field_name) == 2:
                    new_field_names = field_name[1]
                    new_field_names = np.asarray(new_field_names)
                    new_field_names = new_field_names.flatten()
                    for fn in set(new_field_names):
                        fO = frame.fieldOutputs[fn]
                        fO = fO.getSubset(region=inst)
                        n_values = len(fO.values)
                        assert n_values > 0, \
                            ERROR_NO_FIELD_DATA.format(field_name,
                                                       inst_name)

                        fO_location = fO.values[0].position  # NODAL, ELEMENT,
                        if fO_location == NODAL:
                            processPointOutput(fO)
                        elif fO_location in [CENTROID, INTEGRATION_POINT]:
                            cell_data_ = processCellOutput(fO)
                            cell_data = __merge_cellData_dicts(cell_data,
                                                               cell_data_)

                    if fO_location == NODAL:
                        new_point_data = __reshape_fieldOutputs(point_data,
                                                                field_name)
                        point_data.update(new_point_data)
                        for fn in set(new_field_names):
                            del point_data[fn]

                    if fO_location in [CENTROID, INTEGRATION_POINT]:
                        for etype, cd_dict in cell_data.items():
                            new_cd_dict = __reshape_fieldOutputs(cd_dict,
                                                                 field_name)
                            cd_dict.update(new_cd_dict)
                            for fn in set(new_field_names):
                                del cd_dict[fn]

        return points, cells, point_data, cell_data

    if str(type(odbObject)) in ["<type 'OdbInstance'>",
                                "<type 'OdbSet'>"]:
        odbInstance = odbObject
        points, cells, point_data, cell_data = convertInstance(odbInstance,
                                                               frame, 0,
                                                               list_of_outputs)

    elif str(odbObject.__class__) in ["<type 'Odb'>",
                                      "<type 'OdbAssembly'>"]:
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

            if inst_name not in ['Assembly', 'ASSEMBLY', 'assembly']:
                points_, cells_ = convertInstance(inst, frame, idx_shift)[:2]
                points = np.vstack((points, points_))
                cells = __merge_numpy_dicts(cells, cells_)
                idx_shift += len(points_)

    return mo.Mesh(points, cells, point_data, cell_data)


def convertMeshioToMDB(mesh, partname='test', modelname='Model-1', **kwargs):
    print("Warning: only the geometry of the mesh can be converted to ABAQUS´"
          + "model database. Information on field data is lost.")
    assert type(mesh) == mo.mesh.Mesh, 'No meshio Mesh instance'
    all_points = mesh.points
    all_cells = mesh.cells
    n_points = len(all_points)

    # convert nodes to abaqus comaptible representation
    nodeCoords = zip(all_points[:, 0], all_points[:, 1], all_points[:, 2])
    nodeLabels = list(range(1, n_points+1))
    nodeData = [nodeLabels, nodeCoords]

    elementData = []
    element_label_shift = 1
    for etype, els in all_cells.items():
        # convert cells to abaqus compatible representation
        abq_element_type = meshio_to_abaqus_type[etype]
        element_labels = range(element_label_shift,
                               len(els)+element_label_shift)

        # abaqus requires list and int data types
        # consider shifting node labels in connectivity
        element_con = [[int(i+1) for i in x] for x in els]
        elementData.append([abq_element_type, element_labels, element_con])
        # continue counting in next iteration
        element_label_shift += len(els)

    model = mdb.models[modelname]
    partnames = model.parts.keys()
    i = 0
    while partname in partnames:
        if i == 0:
            print('Warning: a part named {} is already in model {}.'.format(partname, modelname))
        partname += '_{}'.format(i)
        i+=1

    model.PartFromNodesAndElements(name=partname, dimensionality=THREE_D,
                                   type=DEFORMABLE_BODY, nodes=nodeData,
                                   elements=elementData)


def convertMeshioToODB(mesh, odbname='test',
                       filename='test.odb', **kwargs):

    all_points = mesh.points
    n_points = len(all_points)
    all_elements = mesh.cells
    all_node_data = mesh.point_data
    all_element_data = mesh.cell_data

    # creat a new odb
    odb = Odb(name=odbname,
              analysisTitle='ODB created from Meshio Instance',
              description='ODB created from Mesio Instance',
              path=filename)

    # add section
    sCat = odb.SectionCategory(name='S5',
                               description='Five-Layered Shell')

    # create part
    odb_part = odb.Part(name='part-1', embeddedSpace=THREE_D,
                        type=DEFORMABLE_BODY)

    # get the nodes
    node_labels = range(1, 1+n_points)
    node_list = zip(node_labels, *[all_points[:, i] for i in (0, 1, 2)])

    # add nodes to odb part
    odb_part.addNodes(nodeData=node_list,
                      nodeSetName='nset-1')

    # get element data
    element_label_shift = 1
    element_labels_LU = {}
    for etype, els in all_elements.items():
        abq_element_type = meshio_to_abaqus_type[etype]
        element_labels = range(element_label_shift,
                               len(els)+element_label_shift)
        # consider shifting node labels in connectivity
        element_con = [[int(i+1) for i in x] for x in els]
        odb_part.addElements(labels=element_labels, connectivity=element_con,
                             type=abq_element_type,
                             elementSetName='{}'.format(etype))
        # continue counting in next iteration
        element_label_shift += len(els)
        element_labels_LU[etype] = element_labels

    # instance part
    odb_inst = odb.rootAssembly.Instance(name='part-1-1',
                                         object=odb_part)

    # create step and frame
    if all_node_data or all_element_data:
        dummy_step = odb.Step(name='step-1', domain=TIME, timePeriod=1.,
                              description='first analysis step')
        dummy_frame = dummy_step.Frame(incrementNumber=0, frameValue=0.,
                                       description='1st frame')

    # write node data
    for nd_name, node_data in all_node_data.items():

        shape = node_data.shape
        if shape == (n_points, ):
            field_type = SCALAR
        elif shape == (n_points, 3):
            field_type = VECTOR
        else:
            print('only scalar and vector data is supported atm')
            continue

        newField = dummy_frame.FieldOutput(name='{}'.format(nd_name),
                                           description='{}'.format(nd_name),
                                           type=field_type)
        if field_type == SCALAR:
            newField.addData(position=CENTROID, instance=odb_inst,
                             labels=element_labels_LU[etype],
                             data=[[x] for x in node_data])

        elif field_type == VECTOR:
            newField.setComponentLabels(('1', '2', '3'))
            newField.setValidInvariants((MAGNITUDE, ))
            newField.addData(position=NODAL, instance=odb_inst,
                             labels=node_labels, data=node_data)
        else:
            print('ERROR processing node output {}'.format(nd_name))

    # write element data
    for etype, ed_dict in all_element_data.items():
        print(etype)
        for ed_name, element_data in ed_dict.items():
            print(ed_name)
            n_elements = len(element_data)
            shape = element_data.shape
            if shape == (n_elements, ):
                field_type = SCALAR
            elif shape == (n_elements, 3):
                field_type = VECTOR
            elif shape == (n_elements, 3, 3):
                field_type = TENSOR_3D_FULL
            else:
                print('only scalar, vector and full3d data is supported atm')
                continue
            if ed_name not in dummy_frame.fieldOutputs.keys():
                # create new
                currentField = (dummy_frame.
                                FieldOutput(name='{}'.format(ed_name),
                                            description='{}'.format(ed_name),
                                            type=field_type))
            else:
                currentField = dummy_frame.fieldOutputs['{}'.format(ed_name)]

            # add data to field_output
            if field_type == SCALAR:
                currentField.addData(position=CENTROID, instance=odb_inst,
                                     labels=element_labels_LU[etype],
                                     data=[[x] for x in element_data])
            elif field_type == VECTOR:
                currentField.setComponentLabels(('1', '2', '3'))
                currentField.setValidInvariants((MAGNITUDE, ))
                currentField.addData(position=CENTROID, instance=odb_inst,
                                     labels=element_labels_LU[etype],
                                     data=[x.tolist() for x in element_data])
            elif field_type == TENSOR_3D_FULL:
                # assume symmetry, keep only the relevant components,
                # reshape to abaqus order
                element_data = [x.flatten()[[0, 4, 8, 1, 2, 5]]
                                for x in element_data]
                currentField.setComponentLabels(('11', '22', '33',
                                                 '12', '13', '23'))
                currentField.setValidInvariants((MISES, TRESCA, PRESS,
                                                 INV3, MAX_PRINCIPAL))
                currentField.addData(position=CENTROID, instance=odb_inst,
                                     labels=element_labels_LU[etype],
                                     data=[x.tolist() for x in element_data])
            else:
                print('ERROR processing element output {}'.format(ed_name))


    pathToODB = odb.path
    odb.save()
    odb.close()
    odb = openOdb(pathToODB)

    return odb
