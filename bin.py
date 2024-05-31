'''
    duplicate_begin_nodes_mask = np.isin(stacked_connections[0], unique_indices, invert=True)
    duplicate_end_nodes_mask = np.isin(stacked_connections[1], unique_indices, invert=True)

    # Find the coordinates of the duplicate begin nodes
    duplicate_begin_nodes_coordinates = rounded_stacked_coordinates[:, stacked_connections[0][duplicate_begin_nodes_mask]]
    duplicate_end_nodes_coordinates = rounded_stacked_coordinates[:, stacked_connections[1][duplicate_end_nodes_mask]]


    print(np.round(unique_nodes.T[0,:], decimals=2))
    print(duplicate_end_nodes_coordinates[0,:])
    indices = np.where(unique_nodes.T[0, :] == duplicate_begin_nodes_coordinates[0, :])[0]

    x_overlap_begin = np.where(np.isin(np.round(unique_nodes.T[0,:], decimals=2), duplicate_begin_nodes_coordinates[0,:]))[0]
    y_overlap_begin = np.where(np.isin(np.round(unique_nodes.T[1, :], decimals=2), duplicate_begin_nodes_coordinates[1, :]))[0]
    z_overlap_begin = np.where(np.isin(np.round(unique_nodes.T[2, :], decimals=2), duplicate_begin_nodes_coordinates[2, :]))[0]

    x_overlap_end = np.where(np.isin(np.round(unique_nodes.T[0, :], decimals=2), duplicate_end_nodes_coordinates[0, :]))[0]
    y_overlap_end = np.where(np.isin(np.round(unique_nodes.T[1, :], decimals=2), duplicate_end_nodes_coordinates[1, :]))[0]
    z_overlap_end = np.where(np.isin(np.round(unique_nodes.T[2, :], decimals=2), duplicate_end_nodes_coordinates[2, :]))[0]



    #a = np.round(unique_nodes.T[0,:], decimals=2)
    #b = duplicate_begin_nodes_coordinates[0,:]
    #print(a)
    #print(b)
    #print()
    #print(np.where(np.isin(a, b))[0])
    #print(np.where(np.all(a == b)))
    #print(x_overlap)
    #print(y_overlap)
    #print(z_overlap)
    print(x_overlap_begin)
    print(y_overlap_begin)
    print(z_overlap_begin)
    print(np.intersect1d(np.intersect1d(x_overlap_begin, y_overlap_begin), z_overlap_begin))

    unique_indices_to_replace_dupl_begins = np.intersect1d(np.intersect1d(x_overlap_begin, y_overlap_begin), z_overlap_begin)
    unique_indices_to_replace_dupl_ends = np.intersect1d(np.intersect1d(x_overlap_end, y_overlap_end), z_overlap_end)

    print(unique_indices_to_replace_dupl_begins.shape, stacked_connections[0][duplicate_begin_nodes_mask].shape)
    stacked_connections[0,:][duplicate_begin_nodes_mask] = unique_indices_to_replace_dupl_begins


    #print('a', duplicate_begin_nodes_mask)
    #print('b', stacked_connections[0][duplicate_begin_nodes_mask])
    #print('c', stacked_coordinates[:, stacked_connections[0][duplicate_begin_nodes_mask]])
    #print('d', np.where(stacked_coordinates == stacked_coordinates[:, stacked_connections[0][duplicate_begin_nodes_mask]]))
    #print(node_indices)
    print('----')

    #print(stacked_coordinates[:, duplicate_begin_nodes_mask])
    remapped_begin_nodes = np.where(stacked_coordinates[:, duplicate_begin_nodes_mask])


    remapped_begin_nodes = np.where(duplicate_begin_nodes_mask, stacked_connections[0], unique_indices[reverse_indices[stacked_connections[0]]])
    remapped_end_nodes = np.where(duplicate_end_nodes_mask, stacked_connections[1], unique_indices[reverse_indices[stacked_connections[1]]])

    # Stack remapped nodes to form remapped connections
    remapped_connections = np.vstack((remapped_begin_nodes, remapped_end_nodes))

    print(remapped_connections)
    print()
    print(self.n_per_hex)
    print(unique_edges.shape)
    '''


'''
        for idx, uc in enumerate(rounded_unique_nodes):
            x_equals = np.where(rounded_stacked_coordinates[0] == uc[0])
            y_equals = np.where(rounded_stacked_coordinates[1] == uc[1])
            z_equals = np.where(rounded_stacked_coordinates[2] == uc[2])
            common_indices = np.intersect1d(np.intersect1d(x_equals, y_equals), z_equals)
            #print(common_indices, common_indices[0])
            #print(np.where(stacked_connections[1] == common_indices[0]))
            print('====')
            stacked_connections[0][np.where(stacked_connections[0] == common_indices)] = idx
            stacked_connections[1][np.where(stacked_connections[1] == common_indices)] = idx


        #print(np.isin(stacked_connections, unique_indices))
        '''

global_coords = self.all_hex_XYZ[0]
global_indices = np.arange(0, self.n_per_hex)

for i in range(1, self.all_hex_connect.shape[0]):
    full_hex_coords = self.all_hex_XYZ[i]
    # print(full_hex_coords)
    # print(global_coords)
    full_hex_connections = self.all_hex_connect[i]
    full_hex_indices = np.arange(i * self.n_per_hex, (i + 1) * self.n_per_hex)
    boolean = np.isin(full_hex_coords, global_coords).all(0)

    b1 = np.isin(global_coords[0], full_hex_coords[0])
    b2 = np.isin(global_coords[0], full_hex_coords[0])
    b3 = np.isin(global_coords[0], full_hex_coords[0])
    b = np.logical_and(b3, np.logical_and(b1, b2))

    # print(global_coords)
    # print(full_hex_coords)
    # print(np.isin(global_coords[0], full_hex_coords[0]))
    # print(global_indices[np.isin(global_coords, full_hex_coords).all(0)])
    full_hex_indices[boolean] = global_indices[boolean]

    # print()
    # for j in range(full_hex_coords.shape[1]):




from Structure_Defs import Geometry_Definition
from helper_functions import calculate_hexagonal_positions
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=7000)

class Hexagonal_Truss(Geometry_Definition):
    def __init__(self, n_rotors: int = 33, r_per_rotor = 12.5, spacing_factor=1, verbose: bool = True):
        self.spacing_factor = spacing_factor
        self.n_rotors = n_rotors
        self.r_rot = r_per_rotor
        self.r_hex = self.r_rot * spacing_factor
        x = int(np.ceil(np.sqrt(self.n_rotors)))
        hex_positions, self.hex_width, self.hex_height, self.hex_area = calculate_hexagonal_positions(n_rotors, self.r_hex, x)
        self.hex_positions = np.array(hex_positions)

        self.X_single_hex = np.array([0, 12.5, 12.5, 0, -12.5, -12.5, 0, 12.5, 12.5, 0, -12.5, -12.5, 0, 0], dtype=float)
        self.Y_single_hex = np.array([0, 0, 0, 0, 0, 0, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 0, 12.5], dtype=float)
        self.Z_single_hex = np.array([-14.434, -7.217, 7.217, 14.434, 7.217, -7.217, -14.434, -7.217, 7.217, 14.434, 7.217, -7.217, 0, 0], dtype=float)

        self.single_hex_mem_idxs = np.array([[0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 4, 5, 7, 8, 10, 11, 5, 2, 12, 1, 5, 0, 3, 6, 9, 2, 4, 5, 1],
                                   [1, 2, 3, 4, 5, 0, 6, 7, 8, 9, 10, 11, 7, 8, 9, 10, 11, 6, 12, 12, 12, 12, 13, 13, 13, 13, 10, 7, 13, 13, 13, 12, 12, 13, 13, 9, 9, 6, 6]])

        self.n_per_hex = self.X_single_hex.size
        self.all_hex_XYZ, self.all_hex_connect = self.transform_coordinates()

        node_indices = np.arange(self.all_hex_XYZ[:,0].ravel().size, dtype=int)

        stacked_coordinates = np.hstack(self.all_hex_XYZ)
        rounded_stacked_coordinates = np.round(stacked_coordinates, decimals=2)
        stacked_connections = np.hstack(self.all_hex_connect)
        unique_nodes, unique_indices = np.unique(rounded_stacked_coordinates.T, axis=0, return_index=True)
        unique_indices = np.arange(unique_indices.size, dtype=int)
        #rounded_unique_nodes = np.round(unique_nodes, decimals=2)
        _, reverse_indices = np.unique(stacked_coordinates.T, axis=0, return_inverse=True)


        coord_to_index_map = {tuple(coord): idx for coord, idx in zip(unique_nodes, unique_indices)}
        new_connections = stacked_connections.copy()

        for i in range(len(node_indices)):
            coord = tuple(rounded_stacked_coordinates.T[i])
            if coord in coord_to_index_map:
                new_index = coord_to_index_map[coord]
                new_connections[stacked_connections == node_indices[i]] = new_index

        unique_edges = np.unique(new_connections, axis=1)
        unique_edges = np.unique(np.sort(unique_edges, axis=0), axis=1)
        print(unique_edges)
        print(np.unique(np.sort(unique_edges, axis=0), axis=1).shape)


        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
        X, Y, Z = unique_nodes[:,0], unique_nodes[:,1], unique_nodes[:,2]
        ax.scatter(X, Y, Z, color='blue', s=10, marker='x')
        plt.show()

        if verbose:
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
        self.plot_structure(show=True)
        super().__init__()

        plt.plot(unique_nodes[:,0], unique_nodes[:,2], linestyle='', marker= 'x')
        for i in range(n_rotors):
            plt.plot(self.all_hex_XYZ[i,0], self.all_hex_XYZ[i,2])
            plt.plot(self.hex_positions[i][0], self.hex_positions[i][1], marker='o', color='red')
        plt.show()
        self.plot_circles(positions=self.hex_positions, width=self.hex_width, height=self.hex_height, title="Hexagonal Layout")

    def transform_coordinates(self):
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
        """Plot circles with given positions."""
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
        fig = plt.figure(figsize=(10,10))
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

        for xz in self.hex_positions:
            ax.scatter(xz[0], avgy, xz[1], color='red')
        if show:
            plt.show()

    def get_XYZ_coords(self):
        return

    def get_member_indices(self):
        return

    def get_section_indices(self):
        return

    def get_bc_indices(self):
        pass

    def get_bc_constraints(self):
        pass

    def get_load_indices(self):
        pass

    def get_applied_loads(self):
        pass

    def get_material_indices(self):
        pass






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




if __name__ == "__main__":
    truss = Hexagonal_Truss(n_rotors=9, r_per_rotor=12.5)
