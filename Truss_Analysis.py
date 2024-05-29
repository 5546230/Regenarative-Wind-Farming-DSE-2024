import numpy as np
import matplotlib.pyplot as plt


class Material:
    def __init__(self, E=190e9, rho=7800, sig_y=340e6):
        'source: [cite]'
        self.E = E     #[Pa]
        self.rho = rho     #[kg/m^3]
        self.sig_y = sig_y  #[Pa]


class Section:
    def __init__(self, radius, thickness,):
        self.R = radius
        self.t = thickness
        self.A = self.calc_area()

    def calc_area(self):
        return 2*np.pi*self.R*self.t


class Mesh:
    def __init__(self, XYZ_coords: np.array, member_indices: np.array, section_ids: np.array, material_ids: np.array, materials, sections):
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

        #print(self.element_Ts.shape, self.element_ks.shape, self.element_Ks.shape)
        #print(self.element_lengths[0])
        #print(self.element_Ks[1])

    def calc_element_lengths(self):
        start_points = self.XYZ_coords[:, member_indices[0, :]]
        end_points = self.XYZ_coords[:, member_indices[1, :]]

        differences = end_points - start_points
        lengths = np.linalg.norm(differences, axis=0)
        return lengths

    def transfer_matrix(self):
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
        single_element_top = np.array([[1, -1], [-1, 1]])[np.newaxis, :,:]
        all_elements_top = np.repeat(single_element_top, self.N_elems, axis=0)

        ks = np.einsum('ijk, i->ijk', all_elements_top, self.element_Es*self.element_As/self.element_lengths)
        return ks

    def transform_stiffness_to_global(self):
        Ks = np.einsum('ikj, ikl, ilm -> ijm', self.element_Ts, self.element_ks, self.element_Ts)
        return Ks

    def plot_structure(self):
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
        plt.show()


class FEM_Solve:
    def __init__(self, mesh: Mesh, bc_indices: np.array, bc_constraints: np.array):
        self.n_dof = 3
        self.mesh = mesh
        self.bc_indices = bc_indices

        constraints = bc_constraints
        global_constraints = np.zeros(mesh.N_nodes*self.n_dof)
        for i, idx in enumerate(self.bc_indices):
            global_constraints[idx*self.n_dof:(idx+1)*self.n_dof] =  constraints[:,i]

        self.active_dofs = mesh.dof_indices[global_constraints==0]
        self.constr_dofs = mesh.dof_indices[global_constraints==1]

        #print(mesh.element_indices)

    def assemble_global(self):
        mesh = self.mesh
        Ks = mesh.element_Ks

        N_active = self.active_dofs.size
        S = np.zeros((N_active, N_active))

        start_points = mesh.element_indices[0, :]
        end_points = mesh.element_indices[1, :]
        dof_range = np.linspace(0, self.n_dof-1, self.n_dof, dtype=int)

        'each column is a node, each row is a x,y,z dof'
        start_dof_idxs = dof_range[:, np.newaxis] + start_points * self.n_dof
        end_dof_idxs = dof_range[:, np.newaxis] +  end_points*self.n_dof

        for i in range(mesh.N_elems):
            elem_start_dofs = start_dof_idxs[:,i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))
            m = np.isin(elem_dofs, self.active_dofs)
            m1, m2 = np.meshgrid(m,m)
            mask  = np.logical_and(m1==True, m2==True)
            S += Ks[i][mask].reshape(S.shape)
        return S


if __name__ == '__main__':
    XYZ_coords = 12*np.array([[-6, 12, 6, -12, 0],
                           [0, 0, 0, 0, 24],
                           [8, 8, -8, -8, 0]])


    member_indices = np.array([[0, 1, 2, 3],
                               [4, 4, 4, 4]])
    section_indices = np.array([0, 0, 0, 0])
    material_indices = np.array([0, 0, 0, 0])

    bc_indices = [0,1,2,3]
    bc_constraints = np.array([[1,1,1,1],
                               [1,1,1,1],
                               [1,1,1,1]])


    steel = Material()
    material_val = Material(E=10000, rho=6600, sig_y=340e6)
    section_val = Section(radius=1, thickness=0.01)
    section_val.A = 8.4
    material_library = [material_val, steel]
    section_library = [section_val, Section(radius=0.5, thickness=0.005)]


    TEST = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices, material_ids=material_indices, materials=material_library, sections=section_library)
    TEST.plot_structure()
    #print(TEST.element_lengths)

    SOLVER = FEM_Solve(mesh=TEST, bc_indices=bc_indices, bc_constraints=bc_constraints)
    SOLVER.assemble_global()