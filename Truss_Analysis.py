import numpy as np
import matplotlib.pyplot as plt


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

    def plot_structure(self):
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
        print(self.active_dofs)
        #print(mesh.element_indices)

    def assemble_global_stiffness(self):
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

    def assemble_loading_vector(self):
        mesh = self.mesh
        global_loading_vector = np.zeros(mesh.N_nodes * self.n_dof)
        for i, idx in enumerate(self.load_indices):
            global_loading_vector[idx*self.n_dof:(idx+1)*self.n_dof] =  self.applied_loads[:,i]
        return global_loading_vector[self.active_dofs]


if __name__ == '__main__':
    XYZ_coords = 12*np.array([[-6, 12, 6, -12, 0],
                           [0, 0, 0, 0, 24],
                           [8, 8, -8, -8, 0]])
    
    def verif_geom():
            XYZ_coords = np.array([[0, 0, 0, 0, 1, 1, 1, 1],
                                    [0, 1, 0, 1, 0, 1, 0, 1],
                                    [0, 0, 1, 1, 0, 0, 1, 1]])
            
            member_indices = np.array([[0,0,0,0,1,1,2,2,3,4,4,4,5,6],
                                       [1,2,3,4,3,5,3,6,7,5,6,7,7,7]])
            
            section_indices = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            material_indices = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0])

            bc_indices = [0,1,4,5]
            bc_constraints= np.array([[1,1,1,1],
                                      [1,0,1,0],
                                      [1,1,1,1]])
            
            return XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints



    member_indices = np.array([[0, 1, 2, 3],
                               [4, 4, 4, 4]])
    section_indices = np.array([0, 0, 0, 0])
    material_indices = np.array([0, 0, 0, 0])

    bc_indices = [0,1,2,3]
    bc_constraints = np.array([[1,1,1,1],
                               [1,1,1,1],
                               [1,1,1,1]])
    

    load_indices = [4]
    applied_loads = np.array([[0],
                              [-100],
                              [-50]])
    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints = verif_geom()

    steel = Material()
    material_val = Material(E=10000, rho=6600, sig_y=340e6)
    section_val = Section(radius=1, thickness=0.01)
    section_val.A = 8.4
    material_library = [material_val, steel]
    section_library = [section_val, Section(radius=0.5, thickness=0.005)]


    TEST = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices, material_ids=material_indices, materials=material_library, sections=section_library)
    TEST.plot_structure()
    #print(TEST.element_lengths)

    SOLVER = FEM_Solve(mesh=TEST, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads)
    S = SOLVER.assemble_global_stiffness()
    P = SOLVER.assemble_loading_vector()
    d = np.linalg.solve(S, P)
