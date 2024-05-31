import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from Structure_Defs import verif_geom_1, tb_val, verif_geom_3, Verif_1, verif_geom_4, hexagon_geom_25
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
from helper_functions import custom_map
np.set_printoptions(linewidth=7000)
'''
NOTES + ASSUMPTIONS:
-  be careful with scaling factor (fix!)
- finalize material library
- organise validation cases
- DO MORE V&V
'''


class Material:
    def __init__(self, E=190e9, rho=7850, sig_y=340e6):
        '''
        :param E: Modulus of Elasticity [Pa]
        :param rho: Density [kg/m^3]
        :param sig_y: yield stress [Pa]
        '''
        self.E = E
        self.rho = rho
        self.sig_y = sig_y


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


    def calc_element_lengths(self):
        '''
        :return: (N_elems, ) Length of each element
        '''
        start_points = self.XYZ_coords[:, self.element_indices[0, :]]
        end_points = self.XYZ_coords[:, self.element_indices[1, :]]

        differences = end_points - start_points
        lengths = np.linalg.norm(differences, axis=0)
        return lengths

    def transfer_matrix(self):
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

    def element_stiffness(self):
        '''
        :return: (N_elems, 2*n_dof, 2*n_dof) multidimensional array containing each local element stiffness
                                                matrix, indexed by the order of self.element_indices
        '''
        single_element_top = np.array([[1, -1], [-1, 1]])[np.newaxis, :,:]
        all_elements_top = np.repeat(single_element_top, self.N_elems, axis=0)
        ks = np.einsum('ijk, i->ijk', all_elements_top, self.element_Es*self.element_As/self.element_lengths)
        return ks

    def element_lumped_mass(self):
        '''
        :return: lumped (2x2) mass matrix, source: [https://www.ce.memphis.edu/7117/notes/presentations/chapter_16.pdf]
        '''
        m_single = np.eye(2)[np.newaxis, :,:]
        m_all = np.repeat(m_single, self.N_elems, axis=0)
        ms = np.einsum('ijk, i->ijk', m_all, self.element_rhos * self.element_As *self.element_lengths/2)
        return ms

    def transform_stiffness_to_global(self, local_matrix):
        '''
        :return: Transformation: local --> global, using tensor product over T
        '''
        Ks = np.einsum('ikj, ikl, ilm -> ijm', self.element_Ts, local_matrix, self.element_Ts)
        return Ks

    def plot_structure(self, show: bool = True):
        '3d plot of nodes and members'
        fig = plt.figure()
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
            plt.plot(Xs, Ys, Zs, color='k')
        if show:
            plt.show()


class FEM_Solve:
    def __init__(self, mesh: Mesh, bc_indices: np.array, bc_constraints: np.array, load_indices: np.array, applied_loads: np.array):
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

    def get_start_end_indices(self):
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

    def assemble_global_stiffness(self):
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

    def assemble_global_mass(self):
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

    def assemble_loading_vector(self):
        '''
        :return: (n_active_nof) Global loading vector sampled over the active indices
        '''
        mesh = self.mesh
        global_loading_vector = np.zeros(mesh.N_nodes * self.n_dof)
        for i, idx in enumerate(self.load_indices):
            global_loading_vector[idx*self.n_dof:(idx+1)*self.n_dof] =  self.applied_loads[:,i]
        return global_loading_vector[self.active_dofs]

    def assemble_self_loading(self):
        '''
        :return: Assemble self-loading (weight) over the active DOFs
        '''

        'assumes gravity acts in z'
        mesh = self.mesh
        global_SL = np.zeros(mesh.N_nodes*self.n_dof)

        'collapse local mass matrices'
        collapsed = np.einsum('ijj->ij', mesh.element_lumped_ms)

        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()
        for i in range(mesh.N_elems):
            elem_start_dofs = start_dof_idxs[:,i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            mask = np.isin(mesh.dof_indices[2::3], elem_dofs)
            global_SL[2::3][mask] -= collapsed[i] * 9.80665
        return global_SL[np.isin(mesh.dof_indices, self.active_dofs)]

    def solve_system(self, factor = 1, include_self_load: bool = False, plot: bool = True,):
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
        M = self.assemble_global_mass()
        K = self.assemble_global_stiffness()
        Minv = np.linalg.inv(M)
        eigenfreqs = np.sqrt(np.abs(np.linalg.eig(-Minv @ K)[0]))
        return np.sort(eigenfreqs)

    def plot_displacements(self, Xp, Yp, Zp):
        '''
        :param Xp: displaced X coordinates
        :param Yp: displaced Y coordinates
        :param Zp: displaced Z coordinates
        '''
        m = self.mesh
        '#d plot of nodes and members'
        fig = plt.figure()
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

            plt.plot(Xs, Ys, Zs, color='k', linestyle='-')
            plt.plot(Xps, Yps, Zps, color='red', linestyle='-')
        plt.show()

    def plot_stresses(self, Xp, Yp, Zp, sigmas):
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

        #norm = plt.Normalize(vmin=min_stress, vmax=max_stress, clip=True)
        norm = TwoSlopeNorm(vmin=min_stress, vcenter=0, vmax=max_stress)
        custom_cmap = custom_map()
        mapper = cm.ScalarMappable(norm=norm, cmap=custom_cmap)#'coolwarm')

        cbar = fig.colorbar(mapper, ax=ax)
        cbar.ax.set_xlabel(r"$\sigma_x$ [MPa]")
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
            legend_string = f"Stress [{i}]: {part_stresses[i]:.2f} MPa"

            ax.plot(Xps, Yps, Zps, color=color, linewidth=.75, label=legend_string)

        handles, labels = ax.get_legend_handles_labels()
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

        legend_ax.legend(handles, labels, loc='right', title="Legend", fontsize=8, ncol=3,bbox_to_anchor=(0.5, 0.5), bbox_transform=plt.gcf().transFigure)
        plt.show()



if __name__ == '__main__':
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
    SOLVER = FEM_Solve(mesh=MESH, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads)
    S = SOLVER.assemble_global_stiffness()

    'solve'
    d, Q, sigma = SOLVER.solve_system(plot=True, factor=100, include_self_load=False)

    print(f'        omega_f [rad/s] = {SOLVER.get_natural_frequencies()}')
    print(f'      displacement [mm] = {d*1000}')
    print(f'   Internal forces [kN] = {Q/1000}')
    print(f'internal stresses [MPa] = {sigma/1e6}')