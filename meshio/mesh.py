# -*- coding: utf-8 -*-
#
import numpy


class Mesh(object):
    def __init__(
        self,
        points,
        cells,
        point_data=None,
        cell_data=None,
        field_data=None,
        node_sets=None,
        gmsh_periodic=None,
    ):
        self.points = points
        self.cells = cells
        self.point_data = point_data if point_data else {}
        self.cell_data = cell_data if cell_data else {}
        self.field_data = field_data if field_data else {}
        self.node_sets = node_sets if node_sets else {}
        self.gmsh_periodic = gmsh_periodic
        return

    def __repr__(self):
        lines = []
        lines.append("Number of points: {}".format(len(self.points)))
        lines.append("Number of elements:")
        for tpe, elems in self.cells.items():
            lines.append("  {}: {}".format(tpe, len(elems)))

        if self.node_sets:
            lines.append("Node sets: {}".format(", ".join(self.node_sets.keys())))

        if self.point_data:
            lines.append("Point data: {}".format(", ".join(self.point_data.keys())))

        cell_data_keys = set()
        for cell_type in self.cell_data:
            cell_data_keys = cell_data_keys.union(self.cell_data[cell_type].keys())
        if cell_data_keys:
            lines.append("Cell data: {}".format(", ".join(cell_data_keys)))

        return "\n".join(lines)

    def prune(self):
        prune_list = []

        for cell_type in ["vertex", "line", "line3"]:
            if cell_type in self.cells:
                prune_list.append(cell_type)

        self.cells.pop("vertex", None)
        self.cells.pop("line", None)
        self.cells.pop("line3", None)
        if "tetra" in self.cells or "tetra10" in self.cells:
            # remove_lower_order_cells
            for cell_type in ["triangle", "triangle6"]:
                if cell_type in self.cells:
                    prune_list.append(cell_type)

        for cell_type in prune_list:
            self.cells.pop(cell_type, None)
            self.cell_data.pop(cell_type, None)

        print("Pruned cell types: {}".format(", ".join(prune_list)))

        # remove_orphaned_nodes.
        # find which nodes are not mentioned in the cells and remove them
        all_cells_flat = numpy.concatenate(
            [vals for vals in self.cells.values()]
        ).flatten()
        orphaned_nodes = numpy.setdiff1d(numpy.arange(len(self.points)), all_cells_flat)
        self.points = numpy.delete(self.points, orphaned_nodes, axis=0)
        # also adapt the point data
        for key in self.point_data:
            self.point_data[key] = numpy.delete(
                self.point_data[key], orphaned_nodes, axis=0
            )

        # reset GLOBAL_ID
        if "GLOBAL_ID" in self.point_data:
            self.point_data["GLOBAL_ID"] = numpy.arange(1, len(self.points) + 1)

        # We now need to adapt the cells too.
        diff = numpy.zeros(len(all_cells_flat), dtype=all_cells_flat.dtype)
        for orphan in orphaned_nodes:
            diff[numpy.argwhere(all_cells_flat > orphan)] += 1
        all_cells_flat -= diff
        k = 0
        for key in self.cells:
            s = self.cells[key].shape
            n = numpy.prod(s)
            self.cells[key] = all_cells_flat[k : k + n].reshape(s)
            k += n
        return
        
    def transform(self, transformation_matrix):
        # transform node coordinates
        homgenous_points = numpy.c_[self.points, numpy.ones(len(self.points))]
        self.points = numpy.einsum('ij,kj -> ki',
                                    transformation_matrix,
                                    homgenous_points)[:,:3]
        
        # transform node fields
        for field_name, field_values in self.point_data.items():
            if len(field_values.shape) == 1:  # SCALAR
                continue
            elif len(field_values.shape) == 2:  # VECTOR
                # get rotation
                R = transformation_matrix[:3, :3]
                rvf = numpy.einsum('ij, kj -> ki',
                                    R, field_values)
                self.point_data.update({field_name: rvf})
            elif len(field_values.shape) == 3:   # TENSOR
                # get rotation
                R = transformation_matrix[:3, :3]
                rtf = numpy.einsum('ij, kjm, nm -> kin',
                                   R, field_values, R)
                self.point_data.update({field_name: rtf})
                
        # transform cell fields
        for etype, cf_dict in self.cell_data.items():
            for field_name, field_values in cf_dict.items():
                if len(field_values.shape) == 1:  # SCALAR
                    continue
                elif len(field_values.shape) == 2:  # VECTOR
                    # get rotation
                    R = transformation_matrix[:3, :3]
                    rvf = numpy.einsum('ij, kj -> ki',
                                        R, field_values)
                    cf_dict.update({field_name: rvf})
                elif len(field_values.shape) == 3:   # TENSOR
                    # get rotation
                    R = transformation_matrix[:3, :3]
                    rtf = numpy.einsum('ij, kjm, nm -> kin',
                                       R, field_values, R)
                    cf_dict.update({field_name: rtf})         
 
 
    def merge(self, other):
        from copy import copy
        points1 = copy(self.points)
        cells1 = copy(self.cells)
        point_data1 = copy(self.point_data)
        cell_data1 = copy(self.cell_data)
        
        points2 = copy(other.points)
        cells2 = copy(other.cells)
        point_data2 = copy(other.point_data)
        cell_data2 = copy(other.cell_data)
        
        # assert(point_data1.keys() == point_data2.keys())
        
        # merge points
        merged_points = numpy.r_[points1, points2]
        # update labels
        for cell_type, cells in cells2.items():
            cells2[cell_type] = len(points1) + cells 
        
        # merge cells
        merged_cells = {}
        etypes1 = cells1.keys() 
        etypes2 = cells2.keys()
             
        # merge cell data
        for etype1 in etypes1:
            if etype1 in etypes2:
                mcells = numpy.r_[cells1[etype1], cells2[etype1]]
                merged_cells[etype1] = mcells
            else:
                # element types only in self
                merged_cells.update({etype1: cells1[etype1]})
               
        for etype2 in etypes2:
            # element types only in other
            if etype2 not in etypes1:
                merged_cells.update({etype2: cells2[etype2]})
                
         # merge point_data
        merged_point_data = {}
        assert point_data1.keys() == point_data2.keys()
        for field_name in point_data1.keys():
            mpd = np.r_[point_data1[field_name], point_data2[field_name]]
            merged_point_data[field_name] = mpd
                
        return Mesh(merged_points, merged_cells, merged_point_data)
