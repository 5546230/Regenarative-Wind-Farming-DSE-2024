from Truss.Honeycomb_Truss import Hexagonal_Truss
from Truss.Truss_Analysis import FEM_Solve, Mesh, Material, Section, Library
from drag_calc import Drag
import numpy as np


class MRS(FEM_Solve):
    def __init__(self, Trated_per_rotor, truss: Hexagonal_Truss,
                 mesh: Mesh, bc_indices: np.array, bc_constraints: np.array, load_indices: np.array, applied_loads: np.array, g_dir: str = 'z'):
        super().__init__(mesh=mesh, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads, g_dir=g_dir)
        self.T_per_rotor = Trated_per_rotor
        self.truss = truss

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
        drag = Drag()
        diams = self.mesh.elem_Ds
        elem_drags = drag.placeholder(L=self.mesh.element_lengths, D=diams)

        elem_nodal_drag = np.array([1/2,1/2])[np.newaxis, :]
        all_nodal_drags = np.repeat(elem_nodal_drag, self.mesh.N_elems, axis=0)
        drags = np.einsum('ij, i->ij', all_nodal_drags, elem_drags)

        start_dof_idxs, end_dof_idxs = self.get_start_end_indices()
        for i in range(self.mesh.N_elems):
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
        return

    def get_inertial_loading(self):
        '''
        :return: inertial loading vector
        '''
        return

    def get_thrust_loading(self):
        T_rated_per_rot = self.T_per_rotor
        front_T_indices = self.get_midpoint_indices(side='front')
        back_T_indices = self.get_midpoint_indices(side='back')

        T_force = np.array([[0], [T_rated_per_rot], [0]]) / 2
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
        P = self.assemble_loading_vector() + self.get_thrust_loading() + self.get_drag_loading()
            # + self.get_drag_loading() + self.get_thrust_loading() + self.get_inertial_loading() + self.get_AFC_loading())
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




if __name__ == "__main__":
    print('running analysis')

    'material and section definitions'
    steel = Material(sig_y=50e6, rho=8000)
    standard_section = Section(radius=.005, thickness=0.001)

    'create libraries'
    material_library = [steel, steel, steel, steel]
    section_library = [standard_section, standard_section, standard_section, standard_section]

    # hex = sizing_truss(Hexagonal_Truss(n_rotors = 3, r_per_rotor = 40.1079757687/2*1.05, spacing_factor=1, verbose=False, depth=25))
    hex = Hexagonal_Truss(n_rotors=3, r_per_rotor=40.1079757687 / 2 * 1.05, spacing_factor=1, verbose=False, depth=25)
    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = hex.function()

    'initialise mesh'
    MESH = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices,
                material_ids=material_indices, materials=material_library, sections=section_library)

    'initialise solver'
    SOLVER = MRS(mesh=MESH, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices,
                       applied_loads=applied_loads, Trated_per_rotor=119e3, truss=hex)

    SOLVER.solve_system()
    'Initialise optimiser'
