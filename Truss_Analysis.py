import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from Structure_Defs import verif_geom_1, tb_val, verif_geom_3, Verif_1
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patheffects as pe
np.set_printoptions(linewidth=7000)


class Material:
    def __init__(self, E=190e9, rho=7800, sig_y=340e6):
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
    def __init__(self, XYZ_coords: np.array, member_indices: np.array, section_ids: np.array, material_ids: np.array, materials: list, sections: list):
        '''
        :param XYZ_coords: (3 x N_nodes) coordinates of each node on the mesh
        :param member_indices: (2 x N_elems) connectivity defining the members by node index
        :param section_ids: (N_elems) section property index of each element
        :param material_ids: (N_elems) material property index of each element
        :param materials: list of materials objects
        :param sections: list of sections objects
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

        self.element_Ts = self.transfer_matrix()
        self.element_ks = self.element_stiffness()
        self.element_Ks = self.transform_stiffness_to_global()

        self.dof_indices = np.linspace(0, N_dof*self.N_nodes-1, N_dof*self.N_nodes, dtype=int)


    def calc_element_lengths(self):
        '''
        :return: (N_elems, ) Length of each element
        '''
        start_points = self.XYZ_coords[:, member_indices[0, :]]
        end_points = self.XYZ_coords[:, member_indices[1, :]]

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
        start_points = self.XYZ_coords[:, member_indices[0, :]]
        end_points = self.XYZ_coords[:, member_indices[1, :]]
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
        :return: (N_elems, 2*n_dof, 2*n_dof) multi-dimensional array containing each local element stiffness
                                                matrix, indexed by the order of self.element_indices
        '''
        single_element_top = np.array([[1, -1], [-1, 1]])[np.newaxis, :,:]
        all_elements_top = np.repeat(single_element_top, self.N_elems, axis=0)

        ks = np.einsum('ijk, i->ijk', all_elements_top, self.element_Es*self.element_As/self.element_lengths)
        return ks

    def transform_stiffness_to_global(self):
        '''
        :return: Transformation local --> global using tensor product over T
        '''
        Ks = np.einsum('ikj, ikl, ilm -> ijm', self.element_Ts, self.element_ks, self.element_Ts)
        return Ks

    def plot_structure(self, show: bool = True):
        '#d plot of nodes and members'
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        for i in range(self.element_indices.shape[1]):
            member_ends =  self.element_indices[:,i]
            Xs = self.X_coords[member_ends]
            Ys = self.Y_coords[member_ends]
            Zs = self.Z_coords[member_ends]
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
        #print(self.active_dofs)


    def get_start_end_indices(self):
        mesh=self.mesh
        start_points = mesh.element_indices[0, :]
        end_points = mesh.element_indices[1, :]
        dof_range = np.linspace(0, self.n_dof - 1, self.n_dof, dtype=int)

        'each column is a node, each row is a x,y,z dof'
        start_dof_idxs = dof_range[:, np.newaxis] + start_points * self.n_dof
        end_dof_idxs = dof_range[:, np.newaxis] + end_points * self.n_dof

        return start_dof_idxs, end_dof_idxs

    def assemble_global_stiffness(self):
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

    def assemble_loading_vector(self):
        mesh = self.mesh
        global_loading_vector = np.zeros(mesh.N_nodes * self.n_dof)
        for i, idx in enumerate(self.load_indices):
            global_loading_vector[idx*self.n_dof:(idx+1)*self.n_dof] =  self.applied_loads[:,i]
        return global_loading_vector[self.active_dofs]

    def solve_system(self, factor = 1, plot: bool = True,):
        S = self.assemble_global_stiffness()
        P = self.assemble_loading_vector()
        d = np.linalg.solve(S, P)

        m = self.mesh

        stacked_coords = np.column_stack((m.X_coords, m.Y_coords, m.Z_coords))
        flattened_coords = stacked_coords.ravel()
        global_displacements = np.zeros_like(flattened_coords)

        flattened_coords[np.isin(m.dof_indices, self.active_dofs)] += d*factor
        global_displacements[np.isin(m.dof_indices, self.active_dofs)] += d*factor

        reshaped_array = flattened_coords.reshape(-1, 3)
        X, Y, Z = reshaped_array[:, 0], reshaped_array[:, 1], reshaped_array[:, 2]
        element_Qs, element_sigmas = self.get_internal_loading(global_ds=global_displacements)

        if plot:
            self.plot_displacements(X, Y, Z)
            self.plot_stresses(X,Y,Z,element_sigmas)
        return d, element_Qs, element_sigmas

    def get_internal_loading(self, global_ds):
        'compute local end displacements'
        mesh=self.mesh
        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()

        Qs = []
        for i in range(mesh.N_elems):
            elem_start_dofs = start_dof_idxs[:, i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            v = global_ds[np.isin(mesh.dof_indices, elem_dofs)]
            T = mesh.element_Ts[i]
            u = T @ v
            k =mesh.element_ks[i]
            Q = k @ u
            Qs.append(-1*Q[1])
        Qs = np.array(Qs)
        sigmas = Qs / mesh.element_As
        return Qs, sigmas

    def plot_displacements(self, Xp, Yp, Zp):
        m = self.mesh
        '#d plot of nodes and members'
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        for i in range(m.element_indices.shape[1]):
            member_ends = m.element_indices[:, i]
            Xs = m.X_coords[member_ends]
            Ys = m.Y_coords[member_ends]
            Zs = m.Z_coords[member_ends]

            Xps = Xp[member_ends]
            Yps = Yp[member_ends]
            Zps = Zp[member_ends]

            plt.plot(Xs, Ys, Zs, color='k')
            plt.plot(Xps, Yps, Zps, color='red')
        plt.show()

    def plot_stresses(self, Xp, Yp, Zp, sigmas):
        m = self.mesh

        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        part_stresses = np.array(sigmas) * 1e-6
        min_stress = 0.85 * np.min(part_stresses)
        max_stress = 1.25 * np.max(part_stresses)

        norm = plt.Normalize(vmin=min_stress, vmax=max_stress, clip=True)
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.YlOrRd)

        cbar = fig.colorbar(mapper, ax=ax)
        cbar.ax.set_xlabel(r"$\sigma_x$ [MPa]")

        sorted_indices = np.argsort(part_stresses)

        for i in sorted_indices:
            member_ends = m.element_indices[:, i]
            Xps = Xp[member_ends]
            Yps = Yp[member_ends]
            Zps = Zp[member_ends]

            color = mapper.to_rgba(part_stresses[i])
            legend_string = f"Stress [{i}]: {part_stresses[i]:.2f} MPa"

            ax.plot(Xps, Yps, Zps, color=color, linewidth=2,
                    path_effects=[pe.Stroke(linewidth=5, foreground='black'), pe.Normal()], label=legend_string)

        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, bbox_to_anchor=(0.5, -0.05), loc='lower center', ncol=3)

        plt.show()


if __name__ == '__main__':

    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = verif_geom_3()

    steel = Material()
    material_val = Material(E=10000, rho=6600, sig_y=340e6)
    section_val = Section(radius=1, thickness=0.01)
    section_val.A = 8.4

    material_val_3 = Material(E=200e9, rho=6600, sig_y=340e6)
    section_val_3 = Section(radius=1, thickness=0.01)
    section_val_3.A = 4000e-6

    material_library = [material_val, material_val_3]
    section_library = [section_val, section_val_3]


    TEST = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices, material_ids=material_indices, materials=material_library, sections=section_library)
    TEST.plot_structure()
    #print(TEST.element_lengths)

    SOLVER = FEM_Solve(mesh=TEST, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads)
    S = SOLVER.assemble_global_stiffness()

    P = SOLVER.assemble_loading_vector()
    d = np.linalg.solve(S, P)

    d, Q, sigma = SOLVER.solve_system(plot=True, factor=1)
    print(d*1000)
    print(sigma/1e6)
