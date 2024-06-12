from Truss.Honeycomb_Truss import Hexagonal_Truss
from Truss.Square_Truss import Square_Truss
from Truss.Truss_Analysis import FEM_Solve, Mesh, Material, Section, Library
from drag_calc import Drag
from Truss.Optimiser import Optimizer
from Truss.helper_functions import CsvOutput
import numpy as np

def cartesian_product(x, y):
    return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])

class MRShex(FEM_Solve):
    def __init__(self, truss_config: dict,
                 mesh: Mesh, bc_indices: np.array, bc_constraints: np.array, load_indices: np.array, applied_loads: np.array, g_dir: str = 'z'):
        '''
        :param truss_config: dict specifying the MRS and loading configuration
        NOTES:
            - extension of FEM_Solve class to include load cases specific to the MRS
        '''
        super().__init__(mesh=mesh, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads, g_dir=g_dir)

        self.T_per_rotor = truss_config['T_per_rotor'][0]
        self.truss = truss_config['truss']
        self.wing_layer_indices = np.array(truss_config['wing_layer_indices'])
        self.drag_per_wing = np.array(truss_config['wing_drags'])
        self.lift_per_wing = np.array(truss_config['wing_lifts'])
        self.moment_per_wing = np.array(truss_config['wing_moments'])
        self.drag_calc = truss_config['drag_calculator']
        self.M_rna = truss_config['M_RNA']
        self.alpha = truss_config['max_alpha_ccwp']

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
        xs = truss.hex_positions[:, 0]
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

    def get_drag_loading(self, multiplier=1):
        '''
        :return: drag force loading vector
        '''
        global_SL = np.zeros(self.mesh.N_nodes * self.n_dof)
        diams = self.mesh.elem_Ds

        'drags only calculated for frontal plane members (which are equal in length)'
        elem_drags = self.drag_calc.placeholder(d=diams)*multiplier
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
        apl_yL = np.max(mesh.Y_coords)
        apl_yM = np.min(mesh.Y_coords)

        indices = []
        loads = []
        for i, z in enumerate(wing_zs):
            'LIFT and DRAG FORCE'
            y_idsL = self.get_xz_plane_indices(y=apl_yL)
            z_ids = self.get_xy_plane_indices(z=z)
            bool_mask_L = np.isin(y_idsL, z_ids)

            temp_indices = y_idsL[bool_mask_L]
            indices.append(temp_indices)

            current_lift = -1*self.lift_per_wing[i]/temp_indices.size
            current_drag = self.drag_per_wing[i]/temp_indices.size

            current_wing_distr_force = np.array([[0], [current_drag], [current_lift]])
            current_loads_LD = np.repeat(current_wing_distr_force, temp_indices.size, axis=1)
            loads.append(current_loads_LD)

            'AERODYNAMIC MOMENT'
            current_M_force = self.moment_per_wing[i]/(apl_yL-apl_yM)/temp_indices.size
            y_idsM = self.get_xz_plane_indices(y=apl_yM)
            bool_mask_M = np.isin(y_idsM, z_ids)
            temp_indices = y_idsL[bool_mask_M]
            indices.append(temp_indices)

            current_wing_distr_force = np.array([[0], [0], [current_M_force]])
            current_loads_M = np.repeat(current_wing_distr_force, temp_indices.size, axis=1)
            loads.append(current_loads_M)

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

        'z +ve upwards: alpha>0-->ccw, alpha<0-->cw'
        if self.alpha>0:
            pass
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

        element_Qs, element_sigmas, reactions = self.get_internal_loading(global_coords=global_coords)
        if plot:
            self.plot_displacements(X, Y, Z)
            self.plot_stresses(X, Y, Z, element_sigmas / factor)
        return d, element_Qs, element_sigmas, reactions

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



class mrsSquare(FEM_Solve):
    def __init__(self, truss_config: dict,
                 mesh: Mesh, bc_indices: np.array, bc_constraints: np.array, load_indices: np.array, applied_loads: np.array, g_dir: str = 'z'):
        '''
        :param truss_config: dict specifying the MRS and loading configuration
        NOTES:
            - extension of FEM_Solve class to include load cases specific to the MRS
        '''
        super().__init__(mesh=mesh, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads, g_dir=g_dir)

        self.T_per_rotor = truss_config['T_per_rotor'][0]
        self.truss = truss_config['truss']
        self.wing_layer_indices = np.array(truss_config['wing_layer_indices'])
        self.drag_per_wing = np.array(truss_config['wing_drags'])
        self.lift_per_wing = np.array(truss_config['wing_lifts'])
        self.moment_per_wing = np.array(truss_config['wing_moments'])
        self.front_drag_calc = truss_config['front_drag_calculator']
        self.side_drag_calc = truss_config['side_drag_calculator']
        self.M_rna = truss_config['M_RNA']
        self.alpha = truss_config['max_alpha_ccwp']

        self.n_per_col = 13
        self.n_per_row = 17

        self.rotor_loading = self.get_rotor_loading()
        self.afc_loading = self.get_AFC_loading()

    def get_xz_plane_indices(self, y: float = 0, tolerance = 0.01):
        c_indices = np.where(np.abs(self.mesh.Y_coords - y) < tolerance)[0]
        return c_indices

    def get_xy_plane_indices(self, z: float = 0, tolerance = 0.01):
        c_indices = np.where(np.abs(self.mesh.Z_coords - z) < tolerance)[0]
        return c_indices

    def get_yz_plane_indices(self, x: float = 0, tolerance = 0.01):
        c_indices = np.where(np.abs(self.mesh.X_coords - x) < tolerance)[0]
        return c_indices

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

    def get_drag_loading(self, multiplier=1):
        '''
        :return: drag force loading vector
        '''
        global_SL = np.zeros(self.mesh.N_nodes * self.n_dof)
        diams = self.mesh.elem_Ds
        ls = self.mesh.element_lengths
        front_drags, side_drags = None, None

        if self.front_drag_calc is not None:
            front_elem_drags = self.front_drag_calc.placeholder(d=diams, l=ls, type='front')*multiplier
            front_elem_nodal_drag = np.array([1/2,1/2])[np.newaxis, :]
            front_all_nodal_drags = np.repeat(front_elem_nodal_drag, self.mesh.N_elems, axis=0)
            front_drags = np.einsum('ij, i->ij', front_all_nodal_drags, front_elem_drags)

        if self.side_drag_calc is not None:
            side_elem_drags = self.side_drag_calc.placeholder(d=diams, l=ls, type='side') * multiplier
            side_elem_nodal_drag = np.array([1 / 2, 1 / 2])[np.newaxis, :]
            side_all_nodal_drags = np.repeat(side_elem_nodal_drag, self.mesh.N_elems, axis=0)
            side_drags = np.einsum('ij, i->ij', side_all_nodal_drags, side_elem_drags)

        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()
        front_indices = self.get_yz_plane_indices(x=np.max(self.mesh.X_coords))
        back_indices = self.get_yz_plane_indices(x=np.min(self.mesh.X_coords))

        side_indices = self.get_xz_plane_indices(y=np.min(self.mesh.Y_coords))

        y_planes = np.unique(self.mesh.Y_coords)[1:]
        for y in y_planes:
            side_indices = np.concatenate((side_indices, self.get_xz_plane_indices(y=y)))

        for i in range(self.mesh.N_elems):
            if self.front_drag_calc is not None:
                if np.all(np.isin(self.mesh.element_indices[:,i], front_indices)==True) or np.all(np.isin(self.mesh.element_indices[:,i], back_indices)==True):
                    #print(self.mesh.element_indices[:, i])
                    elem_start_dofs = start_dof_idxs[:, i]
                    elem_end_dofs = end_dof_idxs[:, i]
                    elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

                    mask = np.isin(self.mesh.dof_indices[0::3], elem_dofs)
                    global_SL[0::3][mask] -= front_drags[i]


            if self.side_drag_calc is not None:
                if self.mesh.Y_coords[self.mesh.element_indices[:, i][0]]==self.mesh.Y_coords[self.mesh.element_indices[:, i][1]]: #np.all(np.isin(self.mesh.element_indices[:, i], side_indices) == True) and
                    #print(self.mesh.Y_coords[self.mesh.element_indices[:, i][0]], self.mesh.Y_coords[self.mesh.element_indices[:, i][1]], self.mesh.element_indices[:, i])
                    elem_start_dofs = start_dof_idxs[:, i]
                    elem_end_dofs = end_dof_idxs[:, i]
                    elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

                    mask = np.isin(self.mesh.dof_indices[1::3], elem_dofs)
                    global_SL[1::3][mask] += side_drags[i]
        return global_SL[np.isin(self.mesh.dof_indices, self.active_dofs)]

    def get_AFC_loading(self):
        '''
        :return: AFC force loading vector
        '''

        AFC_top_back = np.array(range(12, 221, 13))
        AFC_top_front = np.array(range(233, 442, 13))
        total_top_indices = np.concatenate((AFC_top_back, AFC_top_front))
        top_front_load = np.array([[-self.drag_per_wing[2]/total_top_indices.size],
                                   [0],
                                   [-self.lift_per_wing[2]/total_top_indices.size]])*np.ones((3,AFC_top_front.size))
        top_back_load = top_front_load

        AFC_middle_back = np.array(range(8, 217, 13))
        AFC_middle_front = np.array(range(229, 438, 13))
        total_middle_indices = np.concatenate((AFC_middle_back, AFC_middle_front))
        middle_front_load = np.array([[-self.drag_per_wing[1]/total_middle_indices.size],
                                      [0],
                                      [-self.lift_per_wing[1]/total_middle_indices.size]])*np.ones((3,AFC_top_front.size))
        middle_back_load = middle_front_load

        AFC_bottom_back = np.array(range(4, 221, 13))
        AFC_bottom_front = np.array(range(225, 434, 13))
        total_bottom_indices = np.concatenate((AFC_bottom_back, AFC_bottom_front))
        bottom_front_load = np.array([[-self.drag_per_wing[0]/total_bottom_indices.size],
                                      [0],
                                      [-self.lift_per_wing[0]/total_bottom_indices.size]])*np.ones((3,AFC_top_front.size))
        bottom_back_load = bottom_front_load

        indices = np.concatenate((total_top_indices, total_middle_indices, total_bottom_indices))
        loads = np.hstack((top_back_load, top_front_load, middle_back_load, middle_front_load, bottom_back_load, bottom_front_load))

        lv = self.arbitrary_loading_vector(indices, loads)
        return lv

    def get_inertial_loading(self):
        '''
        :return: inertial loading vector
        '''
        return

    def get_rotor_loading(self):
        '''
        :return: thrust loading due to each rotor
        '''
        W_per_rotor = self.M_rna*9.80665
        nodal_load = np.array([[-self.T_per_rotor/4], [0], [-W_per_rotor/4]])

        'first rotor staggering:'
        first_idxs = np.array([1,14,222,235])
        first_col_skip = np.array([0, 6, 9, 15])
        first_row_skip = 2*np.arange(6)

        idxs = []
        loads = []
        for i, row_skip in enumerate(first_row_skip):
            indices = np.repeat(first_idxs[np.newaxis,:], first_col_skip.size, axis=0)+first_col_skip[:,np.newaxis]*self.n_per_col+row_skip
            indices = np.hstack(indices)
            idxs.append(indices)
            current_loads = np.ones((3, indices.size))*nodal_load
            loads.append(current_loads)
        first_loads = np.hstack(loads)
        first_indices = np.hstack(idxs)

        'second rotor staggering:'
        second_idxs = np.array([41, 54, 262, 275])
        second_col_skip = np.array([0,9, ])
        second_row_skip = 2*np.arange(5)

        idxs = []
        loads = []
        for i, row_skip in enumerate(second_row_skip):
            indices = np.repeat(second_idxs[np.newaxis, :], second_col_skip.size, axis=0) + second_col_skip[:, np.newaxis] * self.n_per_col + row_skip
            indices = np.hstack(indices)

            idxs.append(indices)
            current_loads = np.ones((3, indices.size)) * nodal_load
            loads.append(current_loads)
        second_loads = np.hstack(loads)
        second_indices = np.hstack(idxs)

        indices = np.concatenate((first_indices, second_indices))
        loads = np.hstack((first_loads, second_loads))
        lv = self.arbitrary_loading_vector(indices=indices, loads=loads)
        return lv

    def solve_system(self, factor=1, include_self_load: bool = False, plot: bool = False, ) -> tuple[np.array, ...]:
        '''
        :param factor: scaling factor to visualise displacements
        :param plot: bool, plot results
        :return: [-]
        '''
        S = self.assemble_global_stiffness()
        P = self.assemble_loading_vector()  + self.get_drag_loading() + self.afc_loading + self.rotor_loading

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

        element_Qs, element_sigmas, reactions = self.get_internal_loading(global_coords=global_coords)
        #'''
        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()
        Qs = []
        for i in range(self.mesh.N_elems):
            elem_start_dofs = start_dof_idxs[:, i]
            elem_end_dofs = end_dof_idxs[:, i]
            elem_dofs = np.concatenate((elem_start_dofs, elem_end_dofs))

            v = global_displacements[np.isin(self.mesh.dof_indices, elem_dofs)]
            T = self.mesh.element_Ts[i]
            u = T @ v
            k = self.mesh.element_ks[i]
            Q = k @ u
            Qs.append(Q[1])
        #print('TEST:', Qs)
        print(np.array(Qs)[:50])
        print(element_Qs[:50])
        #print((np.sign(Qs)[:500]/np.sign(element_Qs)[:500]))
        element_Qs = np.array(Qs)
        element_sigmas = element_Qs / m.element_As
        #'''
        if plot:
            self.plot_displacements(X, Y, Z)
            #self.plot_stresses(X, Y, Z, element_sigmas / factor)
        return d, element_Qs, element_sigmas, reactions

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
    '''
    config={'wing_layer_indices': [0,1,],
            'wing_lifts': [.5e6, .5e6, 4e6,  .5e6, .5e6, .5e6],
            'wing_moments': [1e6, 1e6, 1e6, 1e6, 1e6, 1e6],
            'wing_drags': [.2e5, .2e5, .2e5, .2e5, .2e5, .2e5],
            'T_per_rotor': [119e3, 119e3],
            'drag_calculator': Drag(V=35, rho=1.225, D_truss=1),
            'M_RNA': 10310.8,
            'max_alpha_ccwp': 0.,
            }
    '''
    # Cm = 0.486
    config_dlc_1 = {'wing_layer_indices': [4, 8, 12],
              'wing_lifts': [1.37e6, 1.42e6, 3.7e6], #[11e6, 11e6, 11e6],#
              'wing_moments': [0,0,0],
              'wing_drags': [.1*1.37e6, .1*1.42e6, .1*3.7e6], #[.66e6, .66e6, .66e6],#
              'T_per_rotor': [150e3, 150e3],
              'front_drag_calculator': Drag(V=10.59, rho=1.225, D_truss=1),
              'side_drag_calculator': None, #,Drag(V=100, rho=1.225, D_truss=1),
              'M_RNA': 37127.73826,
              'max_alpha_ccwp': 0.,
              }

    config_dlc_2 = {'wing_layer_indices': [4, 8, 12],
              'wing_lifts': [1.37e6, 1.42e6, 3.7e6], #[11e6, 11e6, 11e6],#
              'wing_moments': [0,0,0],
              'wing_drags': [.1*1.37e6, .1*1.42e6, .1*3.7e6], #[.66e6, .66e6, .66e6],#
              'T_per_rotor': [150e3*10.59/25, 150e3*10.59/25],
              'front_drag_calculator': Drag(V=25, rho=1.225, D_truss=1),
              'side_drag_calculator': None, #,Drag(V=100, rho=1.225, D_truss=1),
              'M_RNA': 37127.73826,
              'max_alpha_ccwp': 0.,
              }

    config_dlc_3 = {'wing_layer_indices': [4, 8, 12],
              'wing_lifts': [0.0001, 0.0001, 0.0001],#[1.37e6, 1.42e6, 3.7e6], #
              'wing_moments': [0,0,0],
              'wing_drags': [0, 0, 0],#[.1*1.37e6, .1*1.42e6, .1*3.7e6], #
              'T_per_rotor': [0, 0],
              'front_drag_calculator':  None,
              'side_drag_calculator': Drag(V=66, rho=1.225, D_truss=1),
              'M_RNA': 37127.73826,
              'max_alpha_ccwp': 0.,
              }

    config_dlc_4 = {'wing_layer_indices': [4, 8, 12],
              'wing_lifts': [11e6, 11e6, 11e6],#[1.37e6, 1.42e6, 3.7e6], #
              'wing_moments': [0,0,0],
              'wing_drags': [.66e6, .66e6, .66e6],#[.1*1.37e6, .1*1.42e6, .1*3.7e6], #
              'T_per_rotor': [0, 0],
              'front_drag_calculator': Drag(V=66, rho=1.225, D_truss=1),
              'side_drag_calculator': None,
              'M_RNA': 37127.73826,
              'max_alpha_ccwp': 0.,
              }

    'material and section definitions'
    steel = Material()
    standard_section = Section(radius=.05, thickness=0.001)

    'create libraries'
    material_library = [steel, steel, steel, steel]
    section_library = [standard_section, standard_section, standard_section, standard_section]

    # hex = sizing_truss(Hexagonal_Truss(n_rotors = 3, r_per_rotor = 40.1079757687/2*1.05, spacing_factor=1, verbose=False, depth=25))
    #hex = Hexagonal_Truss(n_rotors=33, r_per_rotor=30.38 * 1.05, spacing_factor=1, verbose=False, depth=35)
    #XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = hex.function()

    sq = Square_Truss(depth=50, verbose=False)
    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = sq.function()

    'initialise mesh'
    MESH = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices,
                material_ids=material_indices, materials=material_library, sections=section_library)

    'initialise solver'
    config_dlc_1['truss'] = sq
    config_dlc_2['truss'] = sq
    config_dlc_3['truss'] = sq
    config_dlc_4['truss'] = sq

    SOLVER = mrsSquare(mesh=MESH, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices,
                       applied_loads=applied_loads, truss_config=config_dlc_3, )

    print(f'Izz = {SOLVER.calc_moment_of_inertia_Z():.3e}')
    'Initialise optimiser'
    file_number = 'dlc3'
    csv_output = CsvOutput(f'fem_results_{file_number}.csv')

    OPT = Optimizer(mesh=MESH, solver=SOLVER, plot_output=True, verbose=True, minimum_D=0.36, r_t=1/120) #r_t = 1/120 ==> minimum thickness = 2 mm
    M_change, D_change, ts, sigs, elements, solver, _ = OPT.run_optimisation(tolerance=1e-5, output=csv_output,)
    print(f'D_max = {np.max(D_change[1])}')
    print(f'Izz = {SOLVER.calc_moment_of_inertia_Z():.3e}')
    print(solver.get_natural_frequencies()[:15])

    print('=========================================================================================================')
    print(f'\nInitial Mass = {M_change[0]/1000:.2f} [t], Final Mass = {M_change[1]/1000:.2f} [t]')