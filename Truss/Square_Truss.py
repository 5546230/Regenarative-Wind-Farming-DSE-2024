from Truss.Structure_Defs import Geometry_Definition
from Truss.helper_functions import calculate_hexagonal_positions
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=7000)


class Square_Truss(Geometry_Definition):
    def __init__(self, n_rotors: int = 36, r_per_rotor = 12.5, depth = 12.5, verbose: bool = True):
        super().__init__()
        self.n_rotors = n_rotors
        self.depth = depth
        self.r_rot = r_per_rotor
        self.r_hex = self.r_rot

        x = 4   #int(np.ceil(np.sqrt(self.n_rotors)))
        hex_positions, self.hex_width, self.hex_height, self.hex_area = calculate_hexagonal_positions(n_rotors, self.r_hex, x)
        self.hex_positions = np.array(hex_positions)


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