import numpy as np

def verif_geom_1():
    XYZ_coords = np.array([[0, 0, 0, 0, 1, 1, 1, 1],
                           [0, 1, 0, 1, 0, 1, 0, 1],
                           [0, 0, 1, 1, 0, 0, 1, 1]])

    member_indices = np.array([[0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 0, 5, 0, 5],
                               [2, 3, 3, 3, 6, 7, 6, 7, 7, 7, 4, 4, 1, 1]])

    section_indices = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    material_indices = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    bc_indices = [0, 1, 4, 5]  # [0,1,4,5]
    bc_constraints = np.array([[1, 1, 1, 1],
                               [1, 1, 1, 1],
                               [1, 1, 1, 1]])

    load_indices = [5]
    applied_loads = np.array([[0],
                              [-100],
                              [-50]])

    return XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads


def tb_val():
    XYZ_coords = 12 * np.array([[-6, 12, 6, -12, 0],
                                [0, 0, 0, 0, 24],
                                [8, 8, -8, -8, 0]])

    member_indices = np.array([[0, 1, 2, 3],
                               [4, 4, 4, 4]])
    section_indices = np.array([0, 0, 0, 0])
    material_indices = np.array([0, 0, 0, 0])

    bc_indices = [0, 1, 2, 3]
    bc_constraints = np.array([[1, 1, 1, 1],
                               [1, 1, 1, 1],
                               [1, 1, 1, 1]])

    load_indices = [4]
    applied_loads = np.array([[0],
                              [-100],
                              [-50]])

    return XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads


def verif_geom_3():
    XYZ_coords = np.array([[-4, 4, 4, -4, -2, 2, 2, -2],
                           [0,  0, 0, 0,  10, 10, 10, 10],
                           [4,  4,-4, -4, 2,  2, -2, -2]])

    member_indices = np.array([[0,1,2,3,0,1,2, 3, 4,5,6,7],
                               [4,5,6,7,5,6,7, 4, 5,6,7,4]])

    section_indices = np.ones(12, dtype=int)
    material_indices = np.ones(12, dtype=int)

    bc_indices = [0, 1, 2, 3]  # [0,1,4,5]
    bc_constraints = np.array([[1, 1, 1, 1],
                               [1, 1, 1, 1],
                               [1, 1, 1, 1]])

    load_indices = [4, 7, 6, 5]
    applied_loads = 1000*np.array([[45, 45, 0, 0],
                              [-90,-90, -90, -90],
                              [0, 0,    0, 0]])

    return XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads
