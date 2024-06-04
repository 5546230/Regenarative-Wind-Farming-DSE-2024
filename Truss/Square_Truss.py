from Truss.Structure_Defs import Geometry_Definition
from Truss.helper_functions import calculate_hexagonal_positions
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=7000)


class Square_Truss(Geometry_Definition):
    def __init__(self, n_rotors: int = 17, r_per_rotor = 30.38, depth = 12.5, verbose: bool = True):
        super().__init__()
        self.n_rotors = n_rotors
        self.depth = depth
        self.r_rot = r_per_rotor
        self.r_hex = self.r_rot

        x = 6   #int(np.ceil(np.sqrt(self.n_rotors)))
        half_hex_positions, self.hex_width, self.hex_height, self.hex_area = calculate_hexagonal_positions(n_rotors, self.r_hex, x, column_alignment='vertical')
        self.half_hex_positions = np.array(half_hex_positions)

        for i in range(n_rotors):
            plt.plot(self.half_hex_positions[i][0], self.half_hex_positions[i][1], marker='o', color='red')
        plt.show()
        self.plot_circles(positions=self.half_hex_positions, width=self.hex_width, height=self.hex_height, title="Hexagonal Layout")

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

    def plot_circles(self, positions, width, height, title)->None:
        fig, ax = plt.subplots()
        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect('equal', 'box')
        for (x, y) in positions:
            circle = plt.Circle((x, y), self.r_hex, facecolor='lightblue', alpha=0.6)
            circle2 = plt.Circle((x, y), self.r_hex / 1, facecolor='red', alpha=0.6)
            ax.add_patch(circle)
            ax.add_patch(circle2)
        ax.set_title(title)
        ax.set_xlabel("Width [m]")
        ax.set_ylabel("Height [m]")
        plt.show()



if __name__ == "__main__":
    truss = Square_Truss()