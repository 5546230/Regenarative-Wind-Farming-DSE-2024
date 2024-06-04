from Truss.Honeycomb_Truss import Hexagonal_Truss
from Truss.Truss_Analysis import FEM_Solve, Mesh, Material, Section, Library
from drag_calc import Drag
from Truss.Optimiser import Optimizer
from Truss.helper_functions import CsvOutput
import numpy as np



class MRS(FEM_Solve):
    def __init__(self, truss_config: dict,
                 mesh: Mesh, bc_indices: np.array, bc_constraints: np.array, load_indices: np.array, applied_loads: np.array, g_dir: str = 'z'):
        '''
        :param truss_config: dict specifying the MRS and loading configuration
        NOTES:
            - extension of FEM_Solve class to include load cases specific to the MRS
        '''
        super().__init__(mesh=mesh, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads, g_dir=g_dir)

        self.T_per_rotor = truss_config['T_rated_per_rotor']
        self.truss = truss_config['truss']
        self.wing_layer_indices = np.array(truss_config['wing_layer_indices'])
        self.drag_per_wing = np.array(truss_config['wing_drags'])
        self.lift_per_wing = np.array(truss_config['wing_lifts'])
        self.drag_calc = truss_config['drag_calculator']
        self.M_rna = truss_config['M_RNA']
        self.alpha = truss_config['max_alpha']

    def get_xz_plane_indices(self, y: float = 0, tolerance = 0.01):
        c_indices = np.where(np.abs(self.mesh.Y_coords - y) < tolerance)[0]
        return c_indices

    def get_xy_plane_indices(self, z: float = 0, tolerance = 0.01):
        c_indices = np.where(np.abs(self.mesh.Z_coords - z) < tolerance)[0]
        return c_indices

    def get_yz_plane_indices(self, x: float = 0, tolerance = 0.01):
        c_indices = np.where(np.abs(self.mesh.X_coords - x) < tolerance)[0]
        return c_indices

    def get_midpoint_indices(self, side: str = 'front', tolerance=0.01):
        '''
        :param tolerance: floating point comparison tolerance
        :return: indices of unique nodes lying at hexagon centers, front side
        '''
        truss = self.truss
        xs = truss.hex_positions[:,0]
        zs = truss.hex_positions[:, 1]
        ys = np.ones_like(xs) * np.min(self.mesh.Y_coords)
        if side == 'back':
            ys += truss.depth
        coordinates = np.vstack((xs, ys, zs))

        coordinate_norms = np.linalg.norm(coordinates, axis=0)
        global_norms = np.linalg.norm(self.mesh.XYZ_coords, axis=0)
        c_indices = np.where(np.abs(global_norms[:, np.newaxis] - coordinate_norms) < tolerance)[0]
        return np.array(c_indices)

    def arbitrary_loading_vector(self, indices, loads)-> np.array:
        '''
        :return: (n_active_nof) Global loading vector sampled over the active indices
        , arbitrary loads
        '''
        mesh = self.mesh
        global_loading_vector = np.zeros(mesh.N_nodes * self.n_dof)
        for i, idx in enumerate(indices):
            global_loading_vector[idx*self.n_dof:(idx+1)*self.n_dof] =  loads[:,i]
        return global_loading_vector[self.active_dofs]

    def get_drag_loading(self):
        '''
        :return: drag force loading vector
        '''
        global_SL = np.zeros(self.mesh.N_nodes * self.n_dof)
        diams = self.mesh.elem_Ds

        'drags only calculated for frontal plane members (which are equal in length)'
        elem_drags = self.drag_calc.placeholder(d=diams)
        elem_nodal_drag = np.array([1/2,1/2])[np.newaxis, :]
        all_nodal_drags = np.repeat(elem_nodal_drag, self.mesh.N_elems, axis=0)
        drags = np.einsum('ij, i->ij', all_nodal_drags, elem_drags)

        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()
        front_indices = self.get_xz_plane_indices(y=np.min(self.mesh.Y_coords))
        for i in range(self.mesh.N_elems):
            if np.any(np.isin(self.mesh.element_indices[:,i], front_indices)==False):
                pass
            elem_start_dofs = start_dof_idxs[:, i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            mask = np.isin(self.mesh.dof_indices[1::3], elem_dofs)
            global_SL[1::3][mask] += drags[i]
        return global_SL[np.isin(self.mesh.dof_indices, self.active_dofs)]

    def get_AFC_loading(self):
        '''
        :return: AFC force loading vector
        '''
        mesh = self.mesh

        min_z = np.min(mesh.Z_coords)
        base_apl_z = mesh.Z_coords[np.where(mesh.Z_coords > min_z)][0]
        r_rot = self.truss.r_rot
        next_z_diff = 1.5*r_rot/ np.cos(np.pi/6)
        wing_zs = base_apl_z +  self.wing_layer_indices * next_z_diff
        apl_y = np.max(mesh.Y_coords)

        indices = []
        loads = []
        for i, z in enumerate(wing_zs):
            y_ids = self.get_xz_plane_indices(y=apl_y)
            z_ids = self.get_xy_plane_indices(z=z)
            bool_mask = np.isin(y_ids, z_ids)

            temp_indices = y_ids[bool_mask]
            indices.append(temp_indices)

            current_lift = -1*self.lift_per_wing[i]/temp_indices.size
            current_drag = self.drag_per_wing[i]/temp_indices.size

            current_wing_distr_force = np.array([[0], [current_drag], [current_lift]])
            current_loads = np.repeat(current_wing_distr_force, temp_indices.size, axis=1)
            loads.append(current_loads)

        indices = np.concatenate(indices)
        loads = np.hstack(loads)
        loading_vector = self.arbitrary_loading_vector(indices, loads)
        return loading_vector

    def get_inertial_loading(self):
        '''
        :return: inertial loading vector
        '''
        return

    def get_rotor_loading(self):
        '''
        :return: thrust loading due to each rotor
        '''
        T_rated_per_rot = self.T_per_rotor
        front_T_indices = self.get_midpoint_indices(side='front')
        back_T_indices = self.get_midpoint_indices(side='back')

        T_force = np.array([[0], [T_rated_per_rot], [-self.M_rna*9.80665]]) / 2
        front_T_forces = np.repeat(T_force, front_T_indices.size, axis=1)
        back_T_forces = np.repeat(T_force, back_T_indices.size, axis=1)

        T_forces = np.hstack((front_T_forces, back_T_forces))
        T_indices = np.concatenate((self.get_midpoint_indices(side='front'), self.get_midpoint_indices(side='back')))
        loading_vector = self.arbitrary_loading_vector(indices=T_indices, loads=T_forces)
        return loading_vector

    def solve_system(self, factor=1, include_self_load: bool = False, plot: bool = True, ) -> tuple[np.array, ...]:
        '''
        :param factor: scaling factor to visualise displacements
        :param plot: bool, plot results
        :return: [-]
        '''
        S = self.assemble_global_stiffness()
        P = self.assemble_loading_vector() + self.get_rotor_loading() + self.get_drag_loading() + self.get_AFC_loading()
        if include_self_load:
            P += self.assemble_self_loading()
        d = np.linalg.solve(S, P)

        m = self.mesh
        stacked_coords = np.column_stack((m.X_coords, m.Y_coords, m.Z_coords))
        flattened_coords = stacked_coords.ravel()
        global_displacements = np.zeros_like(flattened_coords)

        flattened_coords[np.isin(m.dof_indices, self.active_dofs)] += d * factor
        global_displacements[np.isin(m.dof_indices, self.active_dofs)] += d * factor

        reshaped_array = flattened_coords.reshape(-1, 3)
        X, Y, Z = reshaped_array[:, 0], reshaped_array[:, 1], reshaped_array[:, 2]
        global_coords = np.vstack((X, Y, Z))

        element_Qs, element_sigmas = self.get_internal_loading(global_coords=global_coords)

        if plot:
            self.plot_displacements(X, Y, Z)
            self.plot_stresses(X, Y, Z, element_sigmas / factor)
        return d, element_Qs, element_sigmas

    def calc_moment_of_inertia_Z(self):
        mesh = self.mesh
        rotational_axis = np.array((np.average(mesh.X_coords),np.average(mesh.Y_coords)))[:, None]

        collapsed_ms = np.einsum('ijj->ij', mesh.element_lumped_ms)
        node_indices = np.arange(mesh.N_nodes, dtype=int)
        global_point_inertias = np.zeros_like(node_indices, dtype=float)
        Xs, Ys = mesh.X_coords, mesh.Y_coords

        for i in range(mesh.N_elems):
            elem_boundaries = mesh.element_indices[:, i]
            mask = np.isin(node_indices, elem_boundaries)

            coords = np.vstack((Xs[mask], Ys[mask]))
            rs = np.linalg.norm(coords-rotational_axis, axis=0)
            global_point_inertias[mask] += collapsed_ms[i] * rs**2
        Izz = np.sum(global_point_inertias)
        return Izz



if __name__ == "__main__":
    config={'wing_layer_indices': [0,],
            'wing_lifts': [4e6, 4e6, 4e6,  4e6, 4e6, 4e6],
            'wing_drags': [1e5, 1e5, 1e5, 1e5, 1e5, 1e5],
            'T_rated_per_rotor': 119e3,
            'drag_calculator': Drag(V=35, rho=1.225, D_truss=1),
            'M_RNA': 10310.8,
            'max_alpha': 0.,
            }

    'material and section definitions'
    steel = Material()
    standard_section = Section(radius=.05, thickness=0.001)

    'create libraries'
    material_library = [steel, steel, steel, steel]
    section_library = [standard_section, standard_section, standard_section, standard_section]

    # hex = sizing_truss(Hexagonal_Truss(n_rotors = 3, r_per_rotor = 40.1079757687/2*1.05, spacing_factor=1, verbose=False, depth=25))
    hex = Hexagonal_Truss(n_rotors=3, r_per_rotor=39.69 / 2 * 1.05, spacing_factor=1, verbose=False, depth=25)
    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = hex.function()

    'initialise mesh'
    MESH = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices,
                material_ids=material_indices, materials=material_library, sections=section_library)

    'initialise solver'
    config['truss'] = hex

    SOLVER = MRS(mesh=MESH, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices,
                       applied_loads=applied_loads, truss_config=config,)

    print(f'Izz = {SOLVER.calc_moment_of_inertia_Z():.3e}')
    'Initialise optimiser'
    file_number = 2
    csv_output = CsvOutput(f'fem_results_{file_number}.csv')

    OPT = Optimizer(mesh=MESH, solver=SOLVER, plot_output=True, verbose=True, minimum_D=0.24) #r_t = 1/120 ==> minimum thickness = 2 mm
    M_change, D_change, ts, sigs, elements = OPT.run_optimisation(tolerance=1e-5, output=csv_output,)

    print('=========================================================================================================')
    print(f'\nInitial Mass = {M_change[0]/1000:.2f} [t], Final Mass = {M_change[1]/1000:.2f} [t]')
    print(D_change)
    print(ts*1000)