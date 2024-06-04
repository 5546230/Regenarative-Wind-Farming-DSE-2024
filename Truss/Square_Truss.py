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

    def __str__(self):
        return f'Square Truss: N_rotors={self.n_rotors}, r_rotor={self.r_rot:.3f} [m], depth={self.depth:.3f} [m]'

    def get_XYZ_coords(self):
        return

    def get_member_indices(self):
        return

    def get_section_indices(self):
        return

    def get_material_indices(self):
        return

    def get_bc_indices(self):
        return self.find_bottom_indices()

    def get_bc_constraints(self):
        return np.ones((3, self.find_bottom_indices().size))

    def get_load_indices(self):
        return

    def get_applied_loads(self):
        return

    def find_bottom_indices(self, tolerance=0.01):
        '''
        :param tolerance: floating point comparison tolerance
        :return: indices of unique nodes lying at hexagon centers, front side
        '''

        z = np.min(self.Z_coords)
        c_indices = np.where(np.abs(self.Z_coords - z) < tolerance)[0]
        return np.array(c_indices)

    @staticmethod
    def get_unique_mesh(all_XYZ, all_connect, node_indices):
        '''
        :param all_XYZ: (n_rotors, 3, n_nodes_per_hex): array containing x,y,z coordinates of all the hexagons
        :param all_connect: (n_rotors, 2, n_edges_per_hex): array defining the connectivity of each hex
        :param node_indices: array defining the indexing of the overlapping mesh
        :return:
        '''

        stacked_coordinates = np.hstack(all_XYZ)
        rounded_stacked_coordinates = np.round(stacked_coordinates, decimals=1)
        stacked_connections = np.hstack(all_connect)
        rounded_unique_nodes, unique_indices = np.unique(rounded_stacked_coordinates.T, axis=0, return_index=True)
        unique_nodes = stacked_coordinates[:, unique_indices].T

        unique_indices = np.arange(unique_indices.size, dtype=int)
        _, reverse_indices = np.unique(stacked_coordinates.T, axis=0, return_inverse=True)

        'map from unique indices to unique nodes'
        coord_to_index_map = {tuple(coord): idx for coord, idx in zip(rounded_unique_nodes, unique_indices)}
        're-map the original connections to the unique node indices'
        new_connections = stacked_connections.copy()

        for i in range(len(node_indices)):
            coord = tuple(rounded_stacked_coordinates.T[i])
            if coord in coord_to_index_map:
                new_index = coord_to_index_map[coord]
                new_connections[stacked_connections == node_indices[i]] = new_index

        'remove any duplicate edges'
        unique_edges = np.unique(new_connections, axis=1)
        unique_edges = np.unique(np.sort(unique_edges, axis=0), axis=1)
        return unique_nodes, unique_indices, unique_edges

    def plot_structure(self, show: bool = True)->None:
        '3d plot of nodes and members'
        fig = plt.figure(figsize=(13,9), layout='constrained')
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
        X, Y, Z = self.X_coords, self.Y_coords, self.Z_coords

        for idx in range(len(X)):
            ax.text(X[idx], Y[idx], Z[idx], f'{idx}', size=8, zorder=1, color='black')

        ax.scatter(X, Y, Z, color='k', s=10)
        bc_indices = self.get_bc_indices()
        mask = np.isin(self.unique_indices, bc_indices)
        ax.scatter(X[mask], Y[mask], Z[mask], color='green', marker='<', s=75, label='Constrained BCs')

        for i in range(self.n_unique_edges):
            member_ends = self.total_edge_con[:, i]
            Xs = X[member_ends]
            Ys = Y[member_ends]
            Zs = Z[member_ends]
            plt.plot(Xs, Ys, Zs, color='k', linewidth=.5)

        avgy = np.average(Y)
        avgx = np.average(X)
        diff = np.max([np.max(X)-np.min(X), np.max(Y)-np.min(Y), np.max(Z)-np.min(Z)])
        ax.set_ylim([avgy-diff/2, avgy+diff/2])
        ax.set_xlim([avgx-diff/2, avgx+diff/2])
        ax.set_zlim([0, diff])

        ax.scatter(self.hex_positions[:,0], avgy, self.hex_positions[:,1], color='red', label='Midpoints')
        ax.set_title(r'$N_{rotors}=$'+f'{self.n_rotors} '+r'$N_p = $'+f'{self.n_unique_nodes}'+r' $N_e = $'+f'{self.n_unique_edges}',fontsize=12)
        if show:
            ax.legend()
            plt.show()

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