from Truss.Honeycomb_Truss import Hexagonal_Truss
from Truss.Truss_Analysis import Mesh, FEM_Solve, Material, Section
from Truss.helper_functions import CsvOutput
import numpy as np

'NEW OPTIMISER CLASS: constant diameter, variable thickness'

class Optimizer:
    def __init__(self, mesh: Mesh, solver: FEM_Solve, minimum_D = 0.1, r_t = 1/120, plot_output: bool = False, verbose: bool = False):
        '''
        :param mesh: instance of Mesh class
        :param solver: instance of Solver class
        :param minimum_D: minimum allowed bar diameter (prevent 0-force members from vanishing)
        :param plot_output: plot before and after opt.
        :param r_t: thickness ratio: r_t := t/D
        :param verbose: print progress
        sources:
            [1]: https://wes.copernicus.org/articles/5/1121/2020/
        '''
        self.mesh = mesh
        self.solver = solver
        self.r_t = r_t
        self.plot_output = plot_output
        self.minimum_D = minimum_D
        self.verbose = verbose

    def calc_buckling_diameter(self, l_i: np.array, E_i: np.array, F_bar_i: np.array, gamma_buckling=1.2):
        '''
        :param l_i: bar lengths
        :param E_i: bar E's
        :param F_bar_i: bar internal axial forces
        :param gamma_buckling: buckling safety factor
        :return: sizing diameter buckling
        '''
        D_i_buckling = (8 * np.abs(F_bar_i) * gamma_buckling * l_i ** 2 / (np.pi ** 3 * E_i * self.r_t)) ** (1 / 4)
        return D_i_buckling

    def calc_yield_diameter(self, F_bar_i: np.array, sigy_i: np.array, gamma_m=1.35):
        '''
        :param F_bar_i: bar internal axial forces
        :param sigy_i: bar yield stresses
        :param gamma_m: material safety factor
        :return: sizing bar diameters for yield
        '''
        D_i = np.sqrt(np.abs(F_bar_i) * gamma_m / (sigy_i * np.pi * self.r_t))
        return D_i

    def calc_A_thin(self, D_i: np.array,):
        '''
        :param D_i: outer diameters
        :return: bar areas (thin walled approx.)
        '''
        A_i = np.pi * D_i ** 2 * self.r_t
        return A_i

    def calc_bar_masses(self, A_i: np.array, l_i: np.array, rho_i: np.array):
        '''
        :param A_i: Cross-sectional areas
        :param l_i: bar lengths
        :param rho_i: bar densities
        :return: bar masses
        '''
        return A_i * l_i * rho_i

    def update_mesh(self, mesh: Mesh, diameters: np.array, )->Mesh:
        '''
        :param mesh: current instance of Mesh
        :param diameters: current array of bar diameters
        :return: mesh with updated geometry
        '''
        mesh.element_As = self.calc_A_thin(D_i=diameters, )
        mesh.element_ks = mesh.element_stiffness()
        mesh.element_lumped_ms = mesh.element_lumped_mass()
        mesh.element_Ks = mesh.transform_stiffness_to_global(local_matrix=mesh.element_ks)
        mesh.element_Ms = mesh.transform_stiffness_to_global(local_matrix=mesh.element_lumped_ms)
        return mesh

    def run_optimisation(self, tolerance = 1e-6, output: CsvOutput = None):
        '''
        :param tolerance: tolerance at which to end diameter opt.
        :param output: (obj) instance of CsvOutput class
        :return: mass (before after), diameters (before after)
        '''

        'initialise mesh'
        m = self.mesh
        s = self.solver
        m.plot_structure()

        'collect element properties'
        Es = m.element_Es
        sigys = np.array([m.materials[i].sig_y for i in m.material_indices])
        Ds = 2 * np.array([m.sections[i].R for i in m.section_indices])
        rhos = m.element_rhos
        lengths = m.element_lengths

        'solve'
        d, Q, sigma = s.solve_system(plot=self.plot_output, factor=1, include_self_load=True)

        D_i_yield = self.calc_yield_diameter(F_bar_i=Q, sigy_i=sigys, )
        D_i_buckling = self.calc_buckling_diameter(l_i=lengths, E_i=Es, F_bar_i=Q, )
        Ds_sizing = np.max(np.vstack((D_i_yield, D_i_buckling)), axis=0)
        Ds_sizing[np.where(Ds_sizing < self.minimum_D)] = self.minimum_D

        if self.verbose:
            print(f'Sigma_max = {np.max(np.abs(sigma) / 1e6)}')
            print('================================================================================')

        D_output = Ds.copy()
        if np.any(Ds < Ds_sizing-tolerance) or np.any(Ds > Ds_sizing+tolerance):
            while np.any(Ds < Ds_sizing-tolerance) or np.any(Ds > Ds_sizing+tolerance):
                if self.verbose:
                    print(f'Sigma_max = {(np.max(np.abs(sigma) / 1e6)):.3f}, max_diff={np.max(Ds_sizing - Ds):.3e},')
                Ds = Ds_sizing
                s.mesh = self.update_mesh(mesh=m, diameters=Ds)
                d, Q, sigma = s.solve_system(plot=False, factor=1, include_self_load=True)

                D_i_yield = self.calc_yield_diameter(F_bar_i=Q, sigy_i=sigys, )
                D_i_buckling = self.calc_buckling_diameter(l_i=lengths, E_i=Es, F_bar_i=Q, )
                Ds_sizing = np.max(np.vstack((D_i_yield, D_i_buckling)), axis=0)
                Ds_sizing[np.where(Ds_sizing < self.minimum_D)] = self.minimum_D

        print(f'\nFinal Sigma_max = {np.max(np.abs(sigma) / 1e6)}')
        D_output = np.vstack((D_output, Ds))
        M_output = [sum(m.elment_total_ms), np.sum(self.calc_bar_masses(A_i=m.element_As, l_i=lengths, rho_i=rhos))]

        Ds = Ds_sizing
        s.mesh =  self.update_mesh(mesh=m, diameters=Ds)

        _,_,sigmas = s.solve_system(plot=self.plot_output, factor=1, include_self_load=True)
        t_output = Ds*self.r_t

        if output is not None:
            print('Writing to output... '+ output.fname)
            fem_output = {
                'Ds [m]': D_output[1],
                'ts [m]': t_output,
                'Sigmas [MPa]': sigmas / 1e6,
                'N_0 [-]': s.mesh.element_indices[0],
                'N_1 [-]': s.mesh.element_indices[1],
            }
            output.write(fem_output)
        return np.array(M_output), D_output, t_output, sigmas, s.mesh.element_indices


def sort_bins(v: np.array, n_bins: int):
    '''
    :param v: array to sort
    :param n_bins: number of bins
    :return: v sorted into n bins
    '''
    min = np.min(v)
    max = np.max(v)

    bins = np.linspace(min, max, n_bins+1)
    bin_placement = np.digitize(v, bins, right=False)
    bin_placement = np.clip(bin_placement - 1, 0, n_bins - 1)
    min_values = bins[1:]

    indices = np.arange(v.size, dtype=int)
    v_sorted = [list(v[bin_placement==i] )for i in range(n_bins)]
    indices_sorted = [list(indices[bin_placement==i] )for i in range(n_bins)]

    return v_sorted, indices_sorted, min_values



if __name__ == "__main__":
    'material and section definitions'
    steel = Material(sig_y=50e6, rho=8000)
    standard_section = Section(radius=.005, thickness=0.001)

    'create libraries'
    material_library = [steel, steel, steel, steel]
    section_library = [standard_section, standard_section, standard_section, standard_section]

    #hex = sizing_truss(Hexagonal_Truss(n_rotors = 3, r_per_rotor = 40.1079757687/2*1.05, spacing_factor=1, verbose=False, depth=25))
    hex = Hexagonal_Truss(n_rotors=3, r_per_rotor=40.1079757687 / 2 * 1.05, spacing_factor=1, verbose=False, depth=25)
    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = hex.function()

    'initialise mesh'
    MESH = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices, material_ids=material_indices, materials=material_library, sections=section_library)

    'initialise solver'
    SOLVER = FEM_Solve(mesh=MESH, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices, applied_loads=applied_loads)

    'Initialise optimiser'
    file_number = 1
    csv_output = CsvOutput(f'fem_results_{file_number}.csv')

    OPT = Optimizer(mesh=MESH, solver=SOLVER, plot_output=True, verbose=True, minimum_D=0.24) #r_t = 1/120 ==> minimum thickness = 2 mm
    M_change, D_change, ts, sigs, elements = OPT.run_optimisation(tolerance=1e-5, output=csv_output,)

    print('=========================================================================================================')
    print(f'\nInitial Mass = {M_change[0]/1000:.2f} [t], Final Mass = {M_change[1]/1000:.2f} [t]')
    print(D_change)
    print(ts*1000)

    #print(sort_bins(v=D_change[1], n_bins=3))