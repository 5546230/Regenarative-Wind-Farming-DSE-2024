from Truss.helper_functions import CsvOutput
import matplotlib.pyplot as plt
from Load_Cases import mrsSquare
from Truss.Truss_Analysis import Mesh, FEM_Solve, Material, Section
from Truss.Square_Truss import Square_Truss
import numpy as np


def mirror_elements(x_start, y_start, z_start, x_end, y_end, z_end, masses, diameters, y_center):
    N = len(masses)

    def mirror_point(x, y, z, y_center):
        return x, 2 * y_center - y, z


    point_dict = {}
    for i in range(N):
        sp = (x_start[i], y_start[i], z_start[i])
        ep = (x_end[i], y_end[i], z_end[i])
        point_dict[(sp, ep)] = i
        point_dict[(ep, sp)] = i

    # Dictionary to keep track of mirrored pairs
    mirrored_pairs = {}

    # Find all mirrored pairs
    for i in range(N):
        xsp_mirrored, ysp_mirrored, zsp_mirrored = mirror_point(x_start[i], y_start[i], z_start[i], y_center)
        xep_mirrored, yep_mirrored, zep_mirrored = mirror_point(x_end[i], y_end[i], z_end[i], y_center)

        sp_mirrored = (xsp_mirrored, ysp_mirrored, zsp_mirrored)
        ep_mirrored = (xep_mirrored, yep_mirrored, zep_mirrored)

        if (sp_mirrored, ep_mirrored) in point_dict:
            j = point_dict[(sp_mirrored, ep_mirrored)]
            mirrored_pairs[i] = j
            mirrored_pairs[j] = i
        elif (ep_mirrored, sp_mirrored) in point_dict:
            j = point_dict[(ep_mirrored, sp_mirrored)]
            mirrored_pairs[i] = j
            mirrored_pairs[j] = i

    # Compare and assign maximum mass to mirrored pairs
    for i, j in mirrored_pairs.items():
        if i < j:  # Ensure each pair is processed only once
            max_mass = max(masses[i], masses[j])
            masses[i] = max_mass
            masses[j] = max_mass

            max_D = max(diameters[i], diameters[j])
            diameters[i] = max_D
            diameters[j] = max_D

    return masses, diameters


def get_new_structure(ms, Ds, Ss, write = False):
    new_ts = Ds / 120
    new_As = np.pi * Ds ** 2 / 120
    steel = Material()
    standard_section = Section(radius=.05, thickness=0.001)
    config_dlc = {'wing_layer_indices': [4, 8, 12],
                  'wing_lifts': [0, 0, 0],  # [1.37e6, 1.42e6, 3.7e6], #
                  'wing_moments': [0, 0, 0],
                  'wing_drags': [0, 0, 0],  # [.1*1.37e6, .1*1.42e6, .1*3.7e6], #
                  'T_per_rotor': [0, 0],
                  'front_drag_calculator': None,
                  'side_drag_calculator': None,
                  'M_RNA': 0,
                  'max_alpha_ccwp': 0.,
                  }

    'create libraries'
    material_library = [steel, steel, steel, steel]
    section_library = [standard_section, standard_section, standard_section, standard_section]

    sq = Square_Truss(depth=50, verbose=False)
    config_dlc['truss'] = sq
    XYZ_coords, member_indices, section_indices, material_indices, bc_indices, bc_constraints, load_indices, applied_loads = sq.function()

    MESH = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, section_ids=section_indices,
                material_ids=material_indices, materials=material_library, sections=section_library)
    MESH.elem_Ds = Ds
    MESH.element_As = new_As
    MESH.element_ks = MESH.element_stiffness()
    MESH.element_lumped_ms = MESH.element_lumped_mass()
    MESH.element_Ks = MESH.transform_stiffness_to_global(local_matrix=MESH.element_ks)
    MESH.element_Ms = MESH.transform_stiffness_to_global(local_matrix=MESH.element_lumped_ms)

    SOLVER = mrsSquare(mesh=MESH, bc_indices=bc_indices, bc_constraints=bc_constraints, load_indices=load_indices,
                       applied_loads=applied_loads, truss_config=config_dlc, )
    print(f'Izz = {SOLVER.calc_moment_of_inertia_Z():.3e}')
    print(SOLVER.get_natural_frequencies()[:15] / (2 * np.pi))

    if write:
        writer = CsvOutput(f'fem_results_final.csv')
        fem_output = {'ms': ms, 'Ds': Ds}
        writer.write(fem_output)
    if Ss is not None:
        plt.plot(Ds, Ss, marker='o', markersize=2, linestyle='', color='k')
        plt.xlabel(r'$D_i$ [m]')
        plt.ylabel(r'$\sigma_i$ [MPa]')
        plt.show()


if __name__ == "__main__":
    reader1 = CsvOutput(f'fem_results_dlc1.csv')
    reader2 = CsvOutput(f'fem_results_dlc2.csv')
    reader3a = CsvOutput(f'fem_results_dlc3.csv')
    reader3b = CsvOutput(f'fem_results_dlc3_right.csv')

    read_data1 = reader1.read()
    read_data2 = reader2.read()
    read_data3a = reader3a.read()
    read_data3b = reader3b.read()

    Ds1 = read_data1['Ds [m]']
    Ss1 = read_data1['Sigmas [MPa]']
    m1 = read_data1['masses [kg]']

    Ds2 = read_data2['Ds [m]']
    Ss2 = read_data2['Sigmas [MPa]']
    m2 = read_data2['masses [kg]']

    Ds3a = read_data3a['Ds [m]']
    Ss3a = read_data3a['Sigmas [MPa]']
    m3a = read_data3a['masses [kg]']

    Ds3b = read_data3b['Ds [m]']
    Ss3b = read_data3b['Sigmas [MPa]']
    m3b= read_data3b['masses [kg]']

    xs_0 = read_data3a['x_0']
    ys_0 = read_data3a['y_0']
    zs_0 = read_data3a['z_0']
    xs_1 = read_data3a['x_1']
    ys_1 = read_data3a['y_1']
    zs_1 = read_data3a['z_1']

    'mirror third case'
    print(f'M3_a, M3_b = {sum(m3a)/1000, sum(m3b)/1000}')
    m3 = np.max(np.vstack((m3a, m3b)), axis=0)
    D3 = np.max(np.vstack((Ds3a, Ds3b)), axis=0)
    print(f'M3_combined = {np.sum(m3)/1000}')
    get_new_structure(ms=m3, Ds=D3, Ss=None, write=False)
    print('-----------------------------------------------------------------------------------------')

    m =  np.max(np.vstack((m3, m1, m2)), axis=0)
    D = np.max(np.vstack((Ds3a, Ds3b, Ds1, Ds2)), axis=0)
    print(f'M_final = {np.sum(m) / 1000}')
    print(f'D_max_final = {np.max(D)}')

    get_new_structure(ms=m, Ds=D, Ss=None, write=True)

