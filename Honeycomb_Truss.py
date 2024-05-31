from Structure_Defs import Geometry_Definition
from helper_functions import calculate_hexagonal_positions
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=7000)


class Hexagonal_Truss(Geometry_Definition):
    def __init__(self, n_rotors: int = 33, r_per_rotor = 12.5, depth = 12.5, spacing_factor=1, verbose: bool = True):
        '''
        :param n_rotors: number of hexagons
        :param r_per_rotor: rotor radius (per MR)
        :param spacing_factor: keep = 1
        :param verbose: print
        '''
        self.spacing_factor = spacing_factor
        self.n_rotors = n_rotors
        self.r_rot = r_per_rotor
        self.r_hex = self.r_rot * spacing_factor
        x = int(np.ceil(np.sqrt(self.n_rotors)))
        hex_positions, self.hex_width, self.hex_height, self.hex_area = calculate_hexagonal_positions(n_rotors, self.r_hex, x)
        self.hex_positions = np.array(hex_positions)

        self.X_single_hex = np.array([0, 12.5, 12.5, 0, -12.5, -12.5, 0, 12.5, 12.5, 0, -12.5, -12.5, 0, 0], dtype=float)*r_per_rotor/12.5
        self.Y_single_hex = np.array([0, 0, 0, 0, 0, 0, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 0, 12.5], dtype=float) * depth/12.5
        self.Z_single_hex = np.array([-14.434, -7.217, 7.217, 14.434, 7.217, -7.217, -14.434, -7.217, 7.217, 14.434, 7.217, -7.217, 0, 0], dtype=float)*r_per_rotor/12.5

        #self.single_hex_mem_idxs = np.array([[0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 4, 5, 7, 8, 10, 11, 5, 2, 12, 1, 5, 0, 3, 6, 9, 2, 4, 5, 1],
        #                           [1, 2, 3, 4, 5, 0, 6, 7, 8, 9, 10, 11, 7, 8, 9, 10, 11, 6, 12, 12, 12, 12, 13, 13, 13, 13, 10, 7, 13, 13, 13, 12, 12, 13, 13, 9, 9, 6, 6]])

        self.single_hex_mem_idxs = np.array([[0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4,   5, 6, 7, 8, 9,  10, 11, 1,  2,  4,  5,  7,  8, 10, 11,  5, 2, 12,  2,  4,  0,  3,  6,  9, 2, 4, 5, 1],
                                             [1, 2, 3, 4, 5, 0, 6, 7, 8, 9, 10, 11, 7, 8, 9, 10, 11, 6, 12, 12, 12, 12, 13, 13, 13, 13, 10, 7, 13, 13, 13, 12, 12, 13, 13, 9, 9, 6, 6]])

        self.n_per_hex = self.X_single_hex.size
        all_hex_XYZ, all_hex_connect = self.transform_coordinates()
        node_indices = np.arange(all_hex_XYZ[:,0].ravel().size, dtype=int)

        stacked_coordinates = np.hstack(all_hex_XYZ)
        rounded_stacked_coordinates = np.round(stacked_coordinates, decimals=2)
        stacked_connections = np.hstack(all_hex_connect)
        unique_nodes, unique_indices = np.unique(rounded_stacked_coordinates.T, axis=0, return_index=True)
        unique_indices = np.arange(unique_indices.size, dtype=int)
        _, reverse_indices = np.unique(stacked_coordinates.T, axis=0, return_inverse=True)

        'map from unique indices to unique nodes'
        coord_to_index_map = {tuple(coord): idx for coord, idx in zip(unique_nodes, unique_indices)}
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

        if verbose:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection='3d')
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')
            ax.set_zlabel('Z [m]')
            X, Y, Z = unique_nodes[:, 0], unique_nodes[:, 1], unique_nodes[:, 2]
            ax.scatter(X, Y, Z, color='blue', s=10, marker='x')
            plt.show()

            print(f'Unique indices: {unique_indices}')
            print(f'number of nodes before = {node_indices.shape[0]}')
            print(f'number of edges before = {stacked_connections.shape[1]}')
            print('----------------------------------------------------------')
            print(f'\nnumber of nodes after = {unique_nodes.shape[0]}')
            print(f'number of edges after = {unique_edges.shape[1]}')

        self.total_edge_con = unique_edges
        self.n_unique_nodes = unique_indices.size
        self.n_unique_edges = unique_edges[0,:].size

        self.X_coords = unique_nodes.T[0,:]
        self.Y_coords = unique_nodes.T[1,:]
        self.Z_coords = unique_nodes.T[2,:]
        super().__init__()

        self.plot_structure(show=True)
        if verbose:
            plt.plot(unique_nodes[:,0], unique_nodes[:,2], linestyle='', marker= 'x')
            for i in range(n_rotors):
                plt.plot(all_hex_XYZ[i,0], all_hex_XYZ[i,2])
                plt.plot(self.hex_positions[i][0], self.hex_positions[i][1], marker='o', color='red')
            plt.show()
            self.plot_circles(positions=self.hex_positions, width=self.hex_width, height=self.hex_height, title="Hexagonal Layout")

    def transform_coordinates(self):
        '''
        :return: re-mapping of single hexagon to a set of all hexagons
        '''
        centre_single_hex_Y = (np.max(self.X_single_hex)+np.min(self.X_single_hex))/2
        centre_single_hex_Z = (np.max(self.X_single_hex)+np.min(self.X_single_hex))/2
        datum_YZ = self.hex_positions[0]

        transformed_X_single = self.X_single_hex + (datum_YZ[0]-centre_single_hex_Y)
        transformed_Y_single = self.Y_single_hex
        transformed_Z_single = self.Z_single_hex + (datum_YZ[1]-centre_single_hex_Z)

        transformation_vectors_xz = self.hex_positions - datum_YZ
        transformed_XYZ_single = np.vstack((transformed_X_single, transformed_Y_single, transformed_Z_single))

        all_hex_coords = np.repeat(transformed_XYZ_single[np.newaxis, :,:], self.n_rotors, axis=0)
        all_hex_coords[:,0] += transformation_vectors_xz[:,0, np.newaxis]
        all_hex_coords[:, 2] += transformation_vectors_xz[:, 1, np.newaxis]

        'transform connectivity'
        all_hex_connectivity = np.repeat(self.single_hex_mem_idxs[np.newaxis, :, :], self.n_rotors, axis=0)
        connectivity_transforms = np.linspace(0, (self.n_rotors-1)*self.n_per_hex, self.n_rotors, dtype=int)

        all_hex_connectivity += connectivity_transforms[:,np.newaxis, np.newaxis]
        return all_hex_coords, all_hex_connectivity


    def plot_circles(self, positions, width, height, title):
        fig, ax = plt.subplots()
        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect('equal', 'box')
        for (x, y) in positions:
            circle = plt.Circle((x, y), self.r_hex, facecolor='lightblue', alpha=0.6)
            circle2 = plt.Circle((x, y), self.r_hex / self.spacing_factor, facecolor='red', alpha=0.6)
            ax.add_patch(circle)
            ax.add_patch(circle2)
        ax.set_title(title)
        ax.set_xlabel("Width [m]")
        ax.set_ylabel("Height [m]")
        plt.show()

    def plot_structure(self, show: bool = True):
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
        ax.set_title(r'$N_{rotors}=$'+f'{self.n_rotors}', fontsize=15)
        for xz in self.hex_positions:
            ax.scatter(xz[0], avgy, xz[1], color='red')
        if show:
            plt.show()

    def get_XYZ_coords(self):
        XYZ_stacked = np.vstack((self.X_coords, self.Y_coords, self.Z_coords))
        return XYZ_stacked

    def get_member_indices(self):
        return self.total_edge_con

    def get_section_indices(self):
        return np.ones(self.n_unique_edges, dtype=int)

    def get_bc_indices(self):
        return np.array([4, 8, 20, 24], dtype=int)

    def get_bc_constraints(self):
        return np.ones((3, 4))

    def get_load_indices(self):
        return np.array([15, 19])

    def get_applied_loads(self):
        return np.array([[0, 1e5],
                         [1e4, 1e4],
                         [0, 0]])

    def get_material_indices(self):
        return np.ones(self.n_unique_edges, dtype=int)

    def function(self):
        return (self.XYZ_array, self.member_indices, self.section_indices, self.material_indices,
                self.bc_indices, self.bc_constraints, self.load_indices, self.applied_loads)



def sizing_truss(n_rotors: int = 33, r_per_rotor = 40.1079757687/2*1.05, depth = 25, spacing_factor=1, verbose: bool = True):
    truss = Hexagonal_Truss(n_rotors, r_per_rotor, depth, spacing_factor, verbose=verbose)

    #truss.load_indices = []
    #truss.applied_loads = []

    truss.bc_indices = np.array([12,56,100,144,188,232,23,67,111,155,109, 243])
    truss.bc_constraints = np.ones((3,truss.bc_indices.shape[1]))

    return truss


if __name__ == "__main__":
    truss = Hexagonal_Truss(n_rotors=33, r_per_rotor=40.1079757687/2*1.05, depth=35)
    #truss = Hexagonal_Truss(n_rotors=1, r_per_rotor=12.5)
    #sizing_truss()