from Truss.Structure_Defs import Geometry_Definition
from Truss.helper_functions import calculate_hexagonal_positions
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=7000)


class Square_Truss(Geometry_Definition):
    def __init__(self):
        super().__init__()

    def get_XYZ_coords(self):
        return

    def get_member_indices(self):
        return

    def get_section_indices(self):
        return

    def get_material_indices(self):
        return

    def get_bc_indices(self):
        return

    def get_bc_constraints(self):
        return

    def get_load_indices(self):
        return

    def get_applied_loads(self):
        return


if __name__ == "__main__":
    truss = Square_Truss()