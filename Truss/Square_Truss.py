from Truss.Structure_Defs import Geometry_Definition
from Truss.helper_functions import calculate_hexagonal_positions
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=7000)


def cartesian_product(x, y):
    return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])


class Square_Truss(Geometry_Definition):
    def __init__(self, n_rotors: int = 17, r_per_rotor = 30.38, depth = 50, n_layers: int = 12, verbose: bool = True):
        self.n_rotors = n_rotors
        self.depth = depth
        self.r_rot = r_per_rotor
        self.r_hex = self.r_rot

        x_cols = 6
        half_hex_positions, self.hex_width, self.hex_height, self.hex_area = calculate_hexagonal_positions(n_rotors, self.r_hex, x_cols, column_alignment='vertical')
        self.half_hex_positions = np.array(half_hex_positions)


        self.X_single_colbox = np.array([depth/2, depth/2, depth/2, depth/2, -depth/2, -depth/2, -depth/2, -depth/2], dtype=float)
        self.Y_single_colbox = np.array([0, 6, 6, 0,  0,  6,  6,  0], dtype=float)#-27.38
        self.Z_single_colbox = np.array([00.00, 00.00, 30.38, 30.38,  00.00,  00.00,  30.38,  30.38], dtype=float)

        self.single_colbox_mem_idxs = np.array([[0, 1, 2, 3, 0,1,2,3,4,4,7,6,0,2,1,4],
                                                [1, 2, 3, 0, 4,5,6,7,5,7,6,5,7,7,6,6]])

        self.X_single_fillbox = np.array([depth/2, depth/2, depth/2, depth/2, -depth/2, -depth/2, -depth/2, -depth/2], dtype=float)
        self.Y_single_fillbox = np.array([0, 23.31, 23.31,0,  0,  23.31,  23.31,  0], dtype=float)
        self.Z_single_fillbox = np.array([00.00, 00.00, 30.38, 30.38,  00.00,  00.00,  30.38,  30.38], dtype=float)

        self.single_fillbox_mem_idxs = np.array([[0, 1, 2, 3, 0,1,2,3,4,4,7,6,0, 1, 5, 2, 3,0,1, 5],
                                                 [1, 2, 3, 0, 4,5,6,7,5,7,6,5,2, 3, 7, 7, 6,5,4, 7]])


        tower_width= 54.76
        self.X_single_tower = np.array([depth / 2, depth / 2, depth / 2, depth / 2, -depth / 2, -depth / 2, -depth / 2, -depth / 2], dtype=float)
        self.Y_single_tower = np.array([0, tower_width/2, tower_width/2, 0, 0, tower_width/2, tower_width/2, 0], dtype=float)
        self.Z_single_tower = np.array([00.00, 00.00, 30.38, 30.38, 00.00, 00.00, 30.38, 30.38], dtype=float)

        self.single_tower_mem_idxs = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 4, 4, 7, 6, 0, 1, 4, 5, 2, 3],
                                                 [1, 2, 3, 0, 4, 5, 6, 7, 5, 7, 6, 5, 2, 3, 6, 7, 7, 6]])

        width_colbox = 6.#abs(np.max(self.Y_single_colbox) - np.min(self.Y_single_colbox))
        height_colbox = 30.38#abs(np.max(self.Z_single_colbox) - np.min(self.Z_single_colbox))
        width_fillbox = 23.31#abs(np.max(self.Y_single_fillbox) - np.min(self.Y_single_fillbox))


        left_width = 3 * width_colbox + 4 * width_fillbox
        print(width_fillbox, width_colbox)

        'LEFT COLBOXES'
        n_left_colboxes = 3*n_layers

        left_colbox_transforms_Y = np.array([0, 2*width_fillbox+width_colbox, 4*width_fillbox+2*width_colbox])
        left_colbox_transforms_Z = np.arange(n_layers, dtype=int) * height_colbox
        left_colbox_trans = cartesian_product(left_colbox_transforms_Y, left_colbox_transforms_Z)

        single_colbox_XYZ = np.vstack((self.X_single_colbox, self.Y_single_colbox, self.Z_single_colbox))
        left_colboxes_XYZ = np.repeat(single_colbox_XYZ[np.newaxis, :,:], n_left_colboxes, axis=0)

        left_colboxes_XYZ[:,1:] += left_colbox_trans[:,:, np.newaxis]
        left_colboxes_connectivity = np.repeat(self.single_colbox_mem_idxs[np.newaxis, :, :], n_left_colboxes, axis=0)
        connectivity_transforms = np.linspace(0, (n_left_colboxes - 1) * self.X_single_colbox.size,  n_left_colboxes, dtype=int)
        left_colboxes_connectivity += connectivity_transforms[:, np.newaxis, np.newaxis]

        node_indices_col = np.arange(left_colboxes_XYZ[:, 0].ravel().size, dtype=int)

        'LEFT FILLBOXES'
        n_left_fillboxes = 4*n_layers

        left_fillbox_transforms_Y = np.array([width_colbox, width_fillbox + width_colbox, 2*width_fillbox + 2*width_colbox, 3*width_fillbox + 2*width_colbox])
        left_fillbox_trans = cartesian_product(left_fillbox_transforms_Y, left_colbox_transforms_Z)
        single_fillbox_XYZ = np.vstack((self.X_single_fillbox, self.Y_single_fillbox, self.Z_single_fillbox))
        left_fillboxes_XYZ = np.repeat(single_fillbox_XYZ[np.newaxis, :, :], n_left_fillboxes, axis=0)

        left_fillboxes_XYZ[:, 1:] += left_fillbox_trans[:, :, np.newaxis]
        left_fillboxes_connectivity = np.repeat(self.single_fillbox_mem_idxs[np.newaxis, :, :], n_left_fillboxes, axis=0)
        connectivity_transforms = np.linspace(0, (n_left_fillboxes - 1) * self.X_single_fillbox.size, n_left_fillboxes, dtype=int)
        left_fillboxes_connectivity += connectivity_transforms[:, np.newaxis, np.newaxis]
        node_indices_fill = np.arange(left_fillboxes_XYZ[:, 0].ravel().size, dtype=int)


        'COMBINING LEFT FILL AND COL BOXES'
        'COL --> FILL'
        node_indices_fill += node_indices_col.size
        node_indices_left = np.concatenate((node_indices_col, node_indices_fill))

        left_fillboxes_connectivity += node_indices_col.size
        left_connectivity = np.hstack((np.hstack(left_colboxes_connectivity), np.hstack(left_fillboxes_connectivity)))

        left_XYZ = np.concatenate((left_colboxes_XYZ, left_fillboxes_XYZ), axis=0)
        left_XYZ = np.hstack(left_XYZ)


        'COMBINING LEFT AND RIGHT AND TOWER'
        node_indices_right = node_indices_left + node_indices_left.size
        right_connectivity = left_connectivity + node_indices_left.size
        right_XYZ = left_XYZ.copy()
        right_XYZ[1,:] += left_width + tower_width


        n_tower = 2*n_layers
        tower_transforms_Y = np.array([left_width, left_width+tower_width/2])
        tower_trans = cartesian_product(tower_transforms_Y, left_colbox_transforms_Z)
        single_tower_XYZ = np.vstack((self.X_single_tower, self.Y_single_tower, self.Z_single_tower))
        tower_XYZ = np.repeat(single_tower_XYZ[np.newaxis, :, :], n_tower, axis=0)
        tower_XYZ[:, 1:] += tower_trans[:, :, np.newaxis]
        tower_connectivity = np.repeat(self.single_tower_mem_idxs[np.newaxis, :, :], n_tower, axis=0)
        connectivity_transforms = np.linspace(0, (n_tower - 1) * self.X_single_tower.size, n_tower,  dtype=int)
        tower_connectivity += connectivity_transforms[:, np.newaxis, np.newaxis]
        node_indices_tower = np.arange(tower_XYZ[:, 0].ravel().size, dtype=int)

        node_indices_tower += node_indices_right.size + node_indices_left.size
        tower_connectivity += node_indices_right.size + node_indices_left.size
        tower_XYZ = np.hstack((tower_XYZ))
        tower_connectivity = np.hstack((tower_connectivity))

        XYZ = np.hstack((left_XYZ, right_XYZ, tower_XYZ))
        node_indices = np.concatenate((node_indices_left, node_indices_right, node_indices_tower))
        connectivity = np.hstack((left_connectivity, right_connectivity, tower_connectivity))

        'Re-map overlapping mesh to unique coordinates'
        unique_nodes, unique_indices, unique_edges = self.get_unique_mesh(XYZ, connectivity, node_indices)


        if verbose:
            print(f'Unique indices: {unique_indices}')
            print(f'number of nodes before = {node_indices_left.shape[0]}')
            print(f'number of edges before = {(left_connectivity).shape[1]}')
            print('----------------------------------------------------------')
            print(f'\nnumber of nodes after = {unique_nodes.shape[0]}')
            print(f'number of edges after = {unique_edges.shape[1]}')


        self.total_edge_con = unique_edges
        self.n_unique_nodes = unique_indices.size
        self.n_unique_edges = unique_edges[0, :].size

        self.X_coords = unique_nodes.T[0, :]
        self.Y_coords = unique_nodes.T[1, :]
        self.Z_coords = unique_nodes.T[2, :]
        self.unique_indices = unique_indices

        self.plot_structure(show=True)
        super().__init__()
        for i in range(n_rotors):
            plt.plot(self.half_hex_positions[i][0], self.half_hex_positions[i][1], marker='o', color='red')
        plt.show()
        self.plot_circles(positions=self.half_hex_positions, width=self.hex_width, height=self.hex_height, title="Hexagonal Layout")

    def __str__(self):
        return f'Square Truss: N_rotors={self.n_rotors*2}, r_rotor={self.r_rot:.3f} [m], depth={self.depth:.3f} [m]'

    def get_XYZ_coords(self):
        XYZ_stacked = np.vstack((self.X_coords, self.Y_coords, self.Z_coords))
        return XYZ_stacked

    def get_member_indices(self):
        return self.total_edge_con

    def get_section_indices(self):
        return np.ones(self.n_unique_edges, dtype=int)

    def get_material_indices(self):
        return np.ones(self.n_unique_edges, dtype=int)

    def get_bc_indices(self):
        return  [312, 325, 338, 91, 104, 117] # self.find_bottom_indices() #

    def get_bc_constraints(self):
        return np.ones((3, self.find_bottom_indices().size))

    def get_load_indices(self):
        return np.array([15, 19, 30, 16, 2])

    def get_applied_loads(self):
        return np.array([[-1e5,-1e5,0, 0,0],
                         [0,0,2e4*0, 2e4*0,2e4*0],
                         [0,0,-3e5*0,-3e5*0,-3e5*0]])

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

        stacked_coordinates = all_XYZ #np.hstack(all_XYZ)
        rounded_stacked_coordinates = np.round(stacked_coordinates, decimals=1)
        stacked_connections = all_connect #np.hstack(all_connect)
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
        fig = plt.figure(figsize=(20,20), layout='constrained')
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

        #ax.scatter(self.hex_positions[:,0], avgy, self.hex_positions[:,1], color='red', label='Midpoints')
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

    def function(self):
        return (self.XYZ_array, self.member_indices, self.section_indices, self.material_indices,
                self.bc_indices, self.bc_constraints, self.load_indices, self.applied_loads)


if __name__ == "__main__":
    truss = Square_Truss(n_layers=12)
    print(truss.find_bottom_indices())