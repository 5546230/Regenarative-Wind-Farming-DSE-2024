from Structure_Defs import Geometry_Definition
import numpy as np
import matplotlib.pyplot as plt

class Hexagonal_Truss(Geometry_Definition):
    def __init__(self, row_distribution = [6, 5, 6, 5, 6, 5]):
        super().__init__()
        self.single_hex_coords = np.array([[0, 12.5, 12.5, 0, -12.5, -12.5, 0, 12.5, 12.5, 0, -12.5, -12.5, 0, 0],
                                                [0, 0, 0, 0, 0, 0, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 0, 12.5],
                                                [-14.434, -7.217, 7.217, 14.434, 7.217, -7.217, -14.434, -7.217, 7.217, 14.434, 7.217, -7.217, 0, 0]],
                                          dtype=float)

        member_indices = np.array([[0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 4, 5, 7, 8, 10, 11, 5,
                                    2, 12, 1, 5, 0, 3, 6, 9, 2, 4, 5, 1],
                                   [1, 2, 3, 4, 5, 0, 6, 7, 8, 9, 10, 11, 7, 8, 9, 10, 11, 6, 12, 12, 12, 12, 13, 13,
                                    13, 13, 10, 7, 13, 13, 13, 12, 12, 13, 13, 9, 9, 6, 6]])

    def visual_check(self):
        return

    def get_XYZ_coords(self):
        return

    def get_member_indices(self):
        return

    def get_section_indices(self):
        return

def hexagon_geom_25():
    '25m diameter hexagon 30m deep'
    XYZ_coords = np.array([[0, 12.5, 12.5, 0, -12.5, -12.5, 0, 12.5, 12.5, 0, -12.5, -12.5, 0, 0],
                           [0, 0, 0, 0, 0, 0, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 0, 12.5],
                           [-14.434, -7.217, 7.217, 14.434, 7.217, -7.217, -14.434, -7.217, 7.217, 14.434, 7.217, -7.217, 0, 0]], dtype=float)

    member_indices = np.array([[0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4,  5,  6, 7, 8, 9, 10, 11,  1,  2,  4,  5,  7,  8, 10, 11,  5, 2, 12,  1,  5,  0,  3,  6,  9, 2, 4, 5, 1],
                               [1, 2, 3, 4, 5, 0, 6, 7, 8, 9, 10, 11, 7, 8, 9, 10, 11, 6, 12, 12, 12, 12, 13, 13, 13, 13, 10, 7, 13, 13, 13, 12, 12, 13, 13, 9, 9, 6, 6]])

    section_indices = np.ones(len(member_indices[0]), dtype=int)
    material_indices = np.ones(len(member_indices[0]), dtype=int)

    bc_indices = [1, 5, 7, 11, 0, 6] #[1, 5, 7, 11, 0, 6]  # [0,1,4,5]
    bc_constraints = np.array([[1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1]])

    load_indices = [12]
    applied_loads = np.array([[0],
                              [-1e5],
                              [-500e2]])

    return XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads