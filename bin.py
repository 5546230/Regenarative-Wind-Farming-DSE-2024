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




from Truss.Structure_Defs import Geometry_Definition
from Truss.helper_functions import calculate_hexagonal_positions
import numpy as np

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




import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from Truss.Structure_Defs import verif_geom_3
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
from Truss.helper_functions import custom_map
from typing import Tuple, Union
np.set_printoptions(linewidth=7000)
'''
NOTES + ASSUMPTIONS:
-  be careful with scaling factor (fix!)
- finalize material library
- organise validation cases
- DO MORE V&V
'''


class Material:
    def __init__(self, E=190e9, rho=7850, sig_y=340e6, sig_el=260e6):
        '''
        :param E: Modulus of Elasticity [Pa]
        :param rho: Density [kg/m^3]
        :param sig_y: yield stress [Pa]
        '''
        self.E = E
        self.rho = rho
        self.sig_y = sig_y
        self.sig_el = sig_el


class Section:
    def __init__(self, radius, thickness,):
        '''
        :param radius: outer radius [m]
        :param thickness: wall thickness [m]
        '''
        self.R = radius
        self.t = thickness
        self.A = self.calc_area()

    def calc_area(self):
        '''
        :return: cs area, thin walled approx. [m^2]
        '''
        return 2*np.pi*self.R*self.t


class Library:
    def __init__(self, elements: Tuple[Union[Material, Section], ...]):
        self.elements = elements

    def add_element(self, element: Union[Material, Section]):
        self.elements += (element,)

    def __iter__(self):
        return iter(self.elements)

    def get_attribute_list(self, attribute: str) -> np.array:
        return np.array([getattr(element, attribute) for element in self.elements])



class Mesh:
    def __init__(self, XYZ_coords: np.array, member_indices: np.array, section_ids: np.array,
                 material_ids: np.array, materials: list, sections: list):
        '''
        :param XYZ_coords: (3 x N_nodes) coordinates of each node on the mesh
        :param member_indices: (2 x N_elems) connectivity defining the members by node index
        :param section_ids: (N_elems) section property index of each element
        :param material_ids: (N_elems) material property index of each element
        :param materials: list of materials objects
        :param sections: list of sections objects
        sources: [https://www.ce.memphis.edu/7117/notes/presentations/chapter_16.pdf]
                [https://repository.bakrie.ac.id/10/1/%5BTSI-LIB-131%5D%5BAslam_Kassimali%5D_Matrix_Analysis_of_Structure.pdf]
        '''

        N_dof = 3
        self.XYZ_coords = XYZ_coords
        self.X_coords = XYZ_coords[0,:]
        self.Y_coords = XYZ_coords[1, :]
        self.Z_coords = XYZ_coords[2, :]

        self.materials = materials
        self.sections = sections

        self.element_indices = member_indices
        self.N_elems = section_ids.size
        self.N_nodes = self.X_coords.size

        self.section_indices = section_ids
        self.material_indices = material_ids
        self.element_lengths = self.calc_element_lengths()

        self.element_Es = np.array([self.materials[i].E for i in self.material_indices])
        self.element_As = np.array([self.sections[i].A for i in self.section_indices])
        self.element_rhos = np.array([self.materials[i].rho for i in self.material_indices])
        self.elment_total_ms = self.element_As * self.element_lengths * self.element_rhos

        self.element_Ts = self.transfer_matrix()
        self.element_ks = self.element_stiffness()
        self.element_lumped_ms = self.element_lumped_mass()
        self.element_Ks = self.transform_stiffness_to_global(local_matrix=self.element_ks)
        self.element_Ms = self.transform_stiffness_to_global(local_matrix=self.element_lumped_ms)

        self.dof_indices = np.linspace(0, N_dof*self.N_nodes-1, N_dof*self.N_nodes, dtype=int)


    def calc_element_lengths(self)->np.array:
        '''
        :return: (N_elems, ) Length of each element
        '''
        start_points = self.XYZ_coords[:, self.element_indices[0, :]]
        end_points = self.XYZ_coords[:, self.element_indices[1, :]]

        differences = end_points - start_points
        lengths = np.linalg.norm(differences, axis=0)
        return lengths

    def transfer_matrix(self)->np.array:
        '''
        :return: (2,n_dof) Transfer matrix T from global to local coordinates
        '''
        I2 = np.eye(2)
        Ts = []
        Ls = self.element_lengths
        start_points = self.XYZ_coords[:, self.element_indices[0, :]]
        end_points = self.XYZ_coords[:, self.element_indices[1, :]]
        differences = end_points - start_points

        cosines = differences / Ls
        for i in range(cosines.shape[1]):
            semi_row = cosines[:,i]
            Ti = np.kron(I2, semi_row)
            Ts.append(Ti)
        Ts = np.array(Ts)
        return Ts

    def element_stiffness(self)->np.array:
        '''
        :return: (N_elems, 2*n_dof, 2*n_dof) multidimensional array containing each local element stiffness
                                                matrix, indexed by the order of self.element_indices
        '''
        single_element_top = np.array([[1, -1], [-1, 1]])[np.newaxis, :,:]
        all_elements_top = np.repeat(single_element_top, self.N_elems, axis=0)
        ks = np.einsum('ijk, i->ijk', all_elements_top, self.element_Es*self.element_As/self.element_lengths)
        return ks

    def element_lumped_mass(self)->np.array:
        '''
        :return: lumped (2x2) mass matrix, source: [https://www.ce.memphis.edu/7117/notes/presentations/chapter_16.pdf]
        '''
        m_single = np.eye(2)[np.newaxis, :,:]
        m_all = np.repeat(m_single, self.N_elems, axis=0)
        ms = np.einsum('ijk, i->ijk', m_all, self.element_rhos * self.element_As *self.element_lengths/2)
        return ms

    def transform_stiffness_to_global(self, local_matrix)->np.array:
        '''
        :return: Transformation: local --> global, using tensor product over T
        '''
        Ks = np.einsum('ikj, ikl, ilm -> ijm', self.element_Ts, local_matrix, self.element_Ts)
        return Ks

    def plot_structure(self, show: bool = True)->None:
        '3d plot of nodes and members'
        fig = plt.figure(figsize=(8,7))
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
        X, Y, Z = self.X_coords, self.Y_coords, self.Z_coords

        for idx in range(len(X)):
            ax.text(X[idx], Y[idx], Z[idx], f'{idx}', size=8, zorder=1, color='black')

        ax.scatter(X,Y,Z, color='k', s=10)
        for i in range(self.element_indices.shape[1]):
            member_ends =  self.element_indices[:,i]
            Xs = X[member_ends]
            Ys = Y[member_ends]
            Zs = Z[member_ends]
            plt.plot(Xs, Ys, Zs, color='k', linewidth=.85)
        if show:
            plt.show()


class FEM_Solve:
    def __init__(self, mesh: Mesh, bc_indices: np.array, bc_constraints: np.array, load_indices: np.array, applied_loads: np.array, g_dir: str = 'z'):
        '''
        :param mesh: Instance of Mesh class defining geometry
        :param bc_indices: node indices on which bcs are applied
        :param bc_constraints: matrix (3, bc_indices.size) specifying the type of constraint on bc nodes
        :param load_indices: node indices on which point loads are applied
        :param applied_loads:
        '''
        self.n_dof = 3
        self.mesh = mesh
        self.bc_indices = bc_indices
        self.load_indices = load_indices
        self.applied_loads = applied_loads

        constraints = bc_constraints
        global_constraints = np.zeros(mesh.N_nodes*self.n_dof)
        for i, idx in enumerate(self.bc_indices):
            global_constraints[idx*self.n_dof:(idx+1)*self.n_dof] =  constraints[:,i]

        self.active_dofs = mesh.dof_indices[global_constraints==0]
        self.constr_dofs = mesh.dof_indices[global_constraints==1]

        dir_dict = {'x':0, 'y': 1, 'z': 2}
        self.dir = dir_dict[g_dir]

    def get_start_end_indices(self) -> tuple[np.array, np.array]:
        '''
        :return: (3 x n_elems) start and end dof indices of each element (3 at start, 3 at end)
        '''
        mesh=self.mesh
        start_points = mesh.element_indices[0, :]
        end_points = mesh.element_indices[1, :]
        dof_range = np.linspace(0, self.n_dof - 1, self.n_dof, dtype=int)

        'each column is a node, each row is a x,y,z dof'
        start_dof_idxs = dof_range[:, np.newaxis] + start_points * self.n_dof
        end_dof_idxs = dof_range[:, np.newaxis] + end_points * self.n_dof

        return start_dof_idxs, end_dof_idxs

    def assemble_global_stiffness(self) -> np.array:
        '''
        :return: (n_active_dof, n_active_dof) global stiffness matrix
        '''
        mesh = self.mesh
        Ks = mesh.element_Ks

        N_active = self.active_dofs.size
        S = np.zeros((N_active, N_active))
        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()

        for i in range(mesh.N_elems):
            elem_start_dofs = start_dof_idxs[:,i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            mask = np.isin(elem_dofs, self.active_dofs)
            S_mask = np.isin(self.active_dofs, elem_dofs)

            S[np.ix_(S_mask, S_mask)] += Ks[i][np.ix_(mask, mask)]
        return S

    def assemble_global_mass(self)-> np.array:
        '''
        :return: (n_active_dof, n_active_dof) global stiffness matrix
        '''
        mesh = self.mesh
        Ms = mesh.element_Ms

        N_active = self.active_dofs.size
        M = np.zeros((N_active, N_active))
        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()

        for i in range(mesh.N_elems):
            elem_start_dofs = start_dof_idxs[:,i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            mask = np.isin(elem_dofs, self.active_dofs)
            M_mask = np.isin(self.active_dofs, elem_dofs)

            M[np.ix_(M_mask, M_mask)] += Ms[i][np.ix_(mask, mask)]
        return M

    def assemble_loading_vector(self)-> np.array:
        '''
        :return: (n_active_nof) Global loading vector sampled over the active indices
        '''
        mesh = self.mesh
        global_loading_vector = np.zeros(mesh.N_nodes * self.n_dof)
        for i, idx in enumerate(self.load_indices):
            global_loading_vector[idx*self.n_dof:(idx+1)*self.n_dof] =  self.applied_loads[:,i]
        return global_loading_vector[self.active_dofs]

    def assemble_self_loading(self)->np.array:
        '''
        :return: Assemble self-loading (weight) over the active DOFs
        '''

        mesh = self.mesh
        global_SL = np.zeros(mesh.N_nodes*self.n_dof)

        'collapse local mass matrices'
        collapsed = np.einsum('ijj->ij', mesh.element_lumped_ms)

        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()
        for i in range(mesh.N_elems):
            elem_start_dofs = start_dof_idxs[:,i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            mask = np.isin(mesh.dof_indices[self.dir::3], elem_dofs)
            global_SL[self.dir::3][mask] -= collapsed[i] * 9.80665
        return global_SL[np.isin(mesh.dof_indices, self.active_dofs)]

    def solve_system(self, factor = 1, include_self_load: bool = False, plot: bool = True,)->tuple[np.array, ...]:
        '''
        :param factor: scaling factor to visualise displacements
        :param plot: bool, plot results
        :return: [-]
        '''
        S = self.assemble_global_stiffness()
        P = self.assemble_loading_vector()
        if include_self_load:
            P += self.assemble_self_loading()
        d = np.linalg.solve(S, P)

        m = self.mesh
        stacked_coords = np.column_stack((m.X_coords, m.Y_coords, m.Z_coords))
        flattened_coords = stacked_coords.ravel()
        global_displacements = np.zeros_like(flattened_coords)

        flattened_coords[np.isin(m.dof_indices, self.active_dofs)] += d*factor
        global_displacements[np.isin(m.dof_indices, self.active_dofs)] += d*factor

        reshaped_array = flattened_coords.reshape(-1, 3)
        X, Y, Z = reshaped_array[:, 0], reshaped_array[:, 1], reshaped_array[:, 2]
        global_coords = np.vstack((X,Y,Z))

        element_Qs, element_sigmas = self.get_internal_loading(global_ds=global_displacements, global_coords=global_coords)

        if plot:
            self.plot_displacements(X, Y, Z)
            self.plot_stresses(X,Y,Z,element_sigmas/factor)
        return d, element_Qs, element_sigmas

    def get_internal_loading(self, global_ds, global_coords):
        '''
        :param global_ds: global displacement vector
        :return: internal forces, internal axial stresses
        '''
        mesh=self.mesh
        #start_dof_idxs, end_dof_idxs = self.get_start_end_indices()

        start_points = global_coords[:, mesh.element_indices[0, :]]
        end_points = global_coords[:, mesh.element_indices[1, :]]

        differences = end_points - start_points
        new_lengths = np.linalg.norm(differences, axis=0)

        delta_length = new_lengths - mesh.element_lengths
        Qs = delta_length * mesh.element_Es * mesh.element_As / mesh.element_lengths
        '''
        #for i in range(mesh.N_elems):
        #    elem_start_dofs = start_dof_idxs[:, i]
        #    elem_end_dofs = end_dof_idxs[:, i]
          #  elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            #v = global_ds[np.isin(mesh.dof_indices, elem_dofs)]
            #T = mesh.element_Ts[i]
            #u = T @ v
            #k =mesh.element_ks[i]
            #Q = k @ u
            #Qs.append(-1*Q[1])

         #   delta_length = new_lengths[i] - mesh.element_lengths[i]
         #   Q_new = delta_length * mesh.element_Es[i]*mesh.element_As[i] / mesh.element_lengths[i]
         #   Qs.append(Q_new)
        '''
        Qs = np.array(Qs)
        sigmas = Qs / mesh.element_As
        return Qs, sigmas

    def get_natural_frequencies(self):
        '''
        :return: array on natural frequencies (n_active_dofs)
        '''
        M = self.assemble_global_mass()
        K = self.assemble_global_stiffness()
        Minv = np.linalg.inv(M)
        eigenfreqs = np.sqrt(np.abs(np.linalg.eig(-Minv @ K)[0]))
        return np.sort(eigenfreqs)

    def plot_displacements(self, Xp, Yp, Zp)->None:
        '''
        :param Xp: displaced X coordinates
        :param Yp: displaced Y coordinates
        :param Zp: displaced Z coordinates
        '''
        m = self.mesh
        '#d plot of nodes and members'
        fig = plt.figure(figsize=(8,7))
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')

        for idx in range(len(Xp)):
            ax.text(Xp[idx], Yp[idx], Zp[idx], f'{idx}', size=8, zorder=1, color='black')

        ax.scatter(Xp, Yp, Zp, color= 'red', s=10)
        for i in range(m.element_indices.shape[1]):
            member_ends = m.element_indices[:, i]
            Xs = m.X_coords[member_ends]
            Ys = m.Y_coords[member_ends]
            Zs = m.Z_coords[member_ends]

            Xps = Xp[member_ends]
            Yps = Yp[member_ends]
            Zps = Zp[member_ends]

            plt.plot(Xs, Ys, Zs, color='k', linestyle='-', linewidth=.85)
            plt.plot(Xps, Yps, Zps, color='red', linestyle='-', linewidth=.85)
        plt.show()

    def plot_stresses(self, Xp, Yp, Zp, sigmas)->None:
        '''
        :param sigmas: axial stresses
        '''
        m = self.mesh

        fig = plt.figure(figsize=(15, 7))
        gs = GridSpec(1, 2, width_ratios=[2.5, 3])
        ax = fig.add_subplot(gs[1], projection='3d')
        legend_ax = fig.add_subplot(gs[0])
        legend_ax.axis('off')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        part_stresses = np.array(sigmas) * 1e-6
        min_stress = 1 * np.min(part_stresses)
        max_stress = 1 * np.max(part_stresses)

        norm = TwoSlopeNorm(vmin=min_stress, vcenter=0, vmax=max_stress)
        custom_cmap = custom_map()
        mapper = cm.ScalarMappable(norm=norm, cmap=custom_cmap)#'coolwarm')

        cbar = fig.colorbar(mapper, ax=ax)
        cbar.ax.set_xlabel(r"$\sigma_\xi$ [MPa]")
        sorted_indices = np.argsort(part_stresses)

        for idx in range(len(Xp)):
            ax.text(Xp[idx], Yp[idx], Zp[idx], f'{idx}', size=12, zorder=100, color='black')

        ax.scatter(Xp, Yp, Zp, color='k')
        for i in sorted_indices:
            member_ends = m.element_indices[:, i]
            Xps = Xp[member_ends]
            Yps = Yp[member_ends]
            Zps = Zp[member_ends]

            color = mapper.to_rgba(part_stresses[i])
            legend_string = f"Stress {i}-({member_ends[0]}, {member_ends[1]}): {part_stresses[i]:.2f} MPa"

            ax.plot(Xps, Yps, Zps, color=color, linewidth=.75, label=legend_string)

        handles, labels = ax.get_legend_handles_labels()
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

        legend_ax.legend(handles, labels, loc='right', title="Legend", fontsize=8, ncol=3,bbox_to_anchor=(0.5, 0.5), bbox_transform=plt.gcf().transFigure)
        plt.show()



if __name__ == '__main__':
    matLib = Library((Material(), Material(rho=6000, E=90e9)))
    secLib = Library((Section(radius=1, thickness=.1), Section(radius=.5, thickness=.01)))
    print(matLib.get_attribute_list('rho'))

    'Geometry definitions'
    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = verif_geom_3()
    #XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = tb_val()
    #XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = verif_geom_1()
    #XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = verif_geom_4()
    #XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = hexagon_geom_25()

    'material and section definitions'
    steel = Material()
    standard_section = Section(radius=0.6, thickness=0.01)

    material_val = Material(E=10000, rho=6600, sig_y=340e6)
    section_val = Section(radius=1, thickness=0.01)
    section_val.A = 8.4

    material_val_3 = Material(E=200e9, rho=6600, sig_y=340e6)
    section_val_3 = Section(radius=.1, thickness=0.0025)
    section_val_3.A = 4000e-6

    material_val_4 = Material(E=70e9, rho=6600, sig_y=340e6)
    section_val_4 = Section(radius=1, thickness=0.01)
    section_val_4.A = 2000e-6

    'create libraries'
    material_library = [material_val, material_val_3, steel, material_val_4]
    section_library = [section_val, section_val_3, standard_section, section_val_4]

    'initialise mesh'
    MESH = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices, material_ids=material_indices, materials=material_library, sections=section_library)
    MESH.plot_structure()

    'initialise solver'
    SOLVER = FEM_Solve(mesh=MESH, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads, g_dir='z')
    S = SOLVER.assemble_global_stiffness()

    'solve'
    d, Q, sigma = SOLVER.solve_system(plot=True, factor=100, include_self_load=False,)

    print(f'        omega_f [rad/s] = {SOLVER.get_natural_frequencies()}')
    print(f'      displacement [mm] = {d*1000}')
    print(f'   Internal forces [kN] = {Q/1000}')
    print(f'internal stresses [MPa] = {sigma/1e6}')


    def get_internal_loading(self, global_coords):
        '''
        :param global_ds: global displacement vector
        :return: internal forces, internal axial stresses
        '''
        mesh = self.mesh
        start_points = global_coords[:, mesh.element_indices[0, :]]
        end_points = global_coords[:, mesh.element_indices[1, :]]

        differences = end_points - start_points
        new_lengths = np.linalg.norm(differences, axis=0)

        delta_length = new_lengths - mesh.element_lengths
        Qs = delta_length * mesh.element_Es * mesh.element_As / mesh.element_lengths
        '''
        #for i in range(mesh.N_elems):
        #    elem_start_dofs = start_dof_idxs[:, i]
        #    elem_end_dofs = end_dof_idxs[:, i]
          #  elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            #v = global_ds[np.isin(mesh.dof_indices, elem_dofs)]
            #T = mesh.element_Ts[i]
            #u = T @ v
            #k =mesh.element_ks[i]
            #Q = k @ u
            #Qs.append(-1*Q[1])

         #   delta_length = new_lengths[i] - mesh.element_lengths[i]
         #   Q_new = delta_length * mesh.element_Es[i]*mesh.element_As[i] / mesh.element_lengths[i]
         #   Qs.append(Q_new)
        '''
        Qs = np.array(Qs)
        sigmas = Qs / mesh.element_As
        return Qs, sigmas






