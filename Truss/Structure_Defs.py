import numpy as np
from abc import ABC, abstractmethod


class Geometry_Definition(ABC):
    def __init__(self):
        self.XYZ_array = self.get_XYZ_coords()
        self.member_indices = self.get_member_indices()
        self.section_indices = self.get_section_indices()
        self.material_indices = self.get_material_indices()
        self.bc_indices = self.get_bc_indices()
        self.bc_constraints = self.get_bc_constraints()
        self.load_indices = self.get_load_indices()
        self.applied_loads = self.get_applied_loads()

    @abstractmethod
    def get_XYZ_coords(self):
        pass

    @abstractmethod
    def get_member_indices(self):
        pass

    @abstractmethod
    def get_section_indices(self):
        pass

    @abstractmethod
    def get_material_indices(self):
        pass

    @abstractmethod
    def get_bc_indices(self):
        pass

    @abstractmethod
    def get_bc_constraints(self):
        pass

    @abstractmethod
    def get_load_indices(self):
        pass

    @abstractmethod
    def get_applied_loads(self):
        pass


class Verif_1(Geometry_Definition):
    def get_XYZ_coords(self):
        return np.array([[0, 0, 0, 0, 1, 1, 1, 1],
                           [0, 1, 0, 1, 0, 1, 0, 1],
                           [0, 0, 1, 1, 0, 0, 1, 1]])

    def get_member_indices(self):
        return np.array([[0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 0, 5, 0, 5],
                               [2, 3, 3, 3, 6, 7, 6, 7, 7, 7, 4, 4, 1, 1]])

    def get_section_indices(self):
        return np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def get_material_indices(self):
        return np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def get_bc_indices(self):
        return np.array([0, 1, 4, 5])

    def get_bc_constraints(self):
        return np.array([[1, 1, 1, 1],
                               [1, 1, 1, 1],
                               [1, 1, 1, 1]])

    def get_load_indices(self):
        return  [5]

    def get_applied_loads(self):
        return np.array([[0],
                              [-100],
                              [-50]])


class Verif_3(Geometry_Definition):
    def get_XYZ_coords(self):
        return np.array([[-4, 4, 4, -4, -2, 2, 2, -2],
                           [0,  0, 0, 0,  10, 10, 10, 10],
                           [4,  4,-4, -4, 2,  2, -2, -2]], dtype=float)

    def get_member_indices(self):
        return np.array([[0,1,2,3,0,1,2, 3, 4,5,6,7],
                               [4,5,6,7,5,6,7, 4, 5,6,7,4]])

    def get_section_indices(self):
        return np.ones(12, dtype=int)

    def get_material_indices(self):
        return np.ones(12, dtype=int)

    def get_bc_indices(self):
        return  [0, 1, 2, 3]

    def get_bc_constraints(self):
        return np.array([[1, 1, 1, 1],
                               [1, 1, 1, 1],
                               [1, 1, 1, 1]])

    def get_load_indices(self):
        return [4, 7, 6, 5]

    def get_applied_loads(self):
        return 1000*np.array([[45, 45, 0, 0],
                              [-90,-90, -90, -90],
                              [0, 0,    0, 0]])


class Verif_Geom_4(Geometry_Definition):
    def get_XYZ_coords(self):
        return np.array([[-3, 9, 2, 0],
                            [0, 0, 0, 12],
                           [5, 6, -3, 0]], dtype=float)

    def get_member_indices(self):
        return np.array([[0, 1, 2],
                               [3, 3, 3]])

    def get_section_indices(self):
        return np.ones(3, dtype=int)*3

    def get_material_indices(self):
        return np.ones(3, dtype=int)*3

    def get_bc_indices(self):
        return   [0, 1, 2]

    def get_bc_constraints(self):
        return np.array([[1, 1, 1],
                               [1, 1, 1],
                               [1, 1, 1]])

    def get_load_indices(self):
        return [3]

    def get_applied_loads(self):
        return 1000*np.array([[75],
                              [-150],
                              [0]])


def verif_geom_1():
    XYZ_coords = np.array([[0, 0, 0, 0, 1, 1, 1, 1],
                           [0, 1, 0, 1, 0, 1, 0, 1],
                           [0, 0, 1, 1, 0, 0, 1, 1]], dtype=float)

    member_indices = np.array([[0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 0, 5, 0, 5, 3, 2, 2],
                               [2, 3, 3, 3, 6, 7, 6, 7, 7, 7, 4, 4, 1, 1, 5, 4, 7]])

    section_indices = np.array(np.ones(len(member_indices[0]), dtype= int))
    material_indices = np.array(np.ones(len(member_indices[0]), dtype=int))

    bc_indices = [0, 1, 4, 5]  # [0,1,4,5]
    bc_constraints = np.array([[1, 1, 1, 1],
                               [1, 0, 1, 0],
                               [1, 1, 1, 1]])

    load_indices = [3, 7]
    applied_loads = np.array([[-1e5, 0],
                              [0, 0],
                              [2e6, 2e6]])

    return XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads

def tb_val():
    'example problem (similar to 8.3)'
    XYZ_coords = 12 * np.array([[-6, 12, 6, -12, 0],
                                [0, 0, 0, 0, 24],
                                [8, 8, -8, -8, 0]], dtype=float)

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
    'problem 8.5'
    XYZ_coords = np.array([[-4, 4, 4, -4, -2, 2, 2, -2],
                           [0,  0, 0, 0,  10, 10, 10, 10],
                           [4,  4,-4, -4, 2,  2, -2, -2]], dtype=float)

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

def verif_geom_4():
    'problem 8.1'
    XYZ_coords = np.array([[-3, 9, 2, 0],
                            [0, 0, 0, 12],
                           [5, 6, -3, 0]], dtype=float)

    member_indices = np.array([[0, 1, 2],
                               [3, 3, 3]])

    section_indices = np.ones(3, dtype=int)*3
    material_indices = np.ones(3, dtype=int)*3

    bc_indices = [0, 1, 2]  # [0,1,4,5]
    bc_constraints = np.array([[1, 1, 1],
                               [1, 1, 1],
                               [1, 1, 1]])

    load_indices = [3]
    applied_loads = 1000*np.array([[75],
                              [-150],
                              [0]])

    return XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads

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