import numpy as np
import matplotlib.pyplot as plt


class OS_Steel():
    def __init__(self):
        self.E = 190e9
        self.rho = 7800
        self.sig_y = 340e6


class Pylon:
    def __init__(self, radius_outer, length, material):
        self.r_out = radius_outer
        self.L = length
        self.mat = material

    def calc_area(self, t):
        '''
        :param t: wall thickness [m]
        :return: CS area
        '''
        return 2 * np.pi * self.r_out * t

    def calc_I(self, t):
        '''
        :param t: thickness [m]
        :return: moment of inertia (thin walled)
        '''
        return np.pi * self.r_out ** 3 * t

    def axial_thickness(self, axial_force):
        '''
        :param axial_force: applied axial force (N)
        :return: wall thickness to sustain axial load within sigma_yield
        '''
        t_a = axial_force / (self.mat.sig_y * 2 * np.pi * self.r_out)
        return t_a

    def bending_thickness(self, point_force, position):
        M = position * point_force
        t_b = M / (np.pi * self.r_out ** 2 * self.mat.sig_y)
        return t_b

    def calc_mass(self, t):
        M = self.calc_area(t=t) * self.L * self.mat.rho
        return M

    def buckling(self, axial_load, n=0.22):
        t = self.L ** 2 * axial_load / (np.pi ** 3 * n * self.mat.E * self.r_out ** 3)
        return t

    def apply_loads(self, load_arr: np.array, pos_arr: np.array, thickness):
        EI = self.mat.E
        EI *= self.calc_I(t=thickness)

        delta_arr = load_arr * pos_arr * pos_arr / (6 * EI) * (3 * self.L - pos_arr)
        return sum(delta_arr)


if __name__ == "__main__":
    R_rot = 29.6  # [m]
    D_rot = 2 * R_rot
    N_rot = 33

    V_infty = 7
    P_rated = 30e6
    T_rated = P_rated / V_infty
    T_perR = T_rated / N_rot

    # W_rna = 38790.46184 * 9.81

    M_truss = 3000e3  # 2235959.595*1.5
    M_RNA = 1280085.241

    print(f'MASS TRUSS: {M_truss / 1000:.4f}\n')

    # ---------------------------------------------------------------------------------
    print('TOWER')

    pyl = Pylon(radius_outer=10, length=60 + 25 + 350, material=OS_Steel())
    M_i = 2900e3
    err = np.infty
    while err > 100:
        t_a = pyl.axial_thickness(axial_force=(M_truss + M_RNA + M_i) * 9.81)
        t_b = pyl.bending_thickness(point_force=T_rated, position=60 + 25 + .5 * 350)
        t_buckl = pyl.buckling(axial_load=(M_truss + M_RNA + M_i) * 9.81)
        M = pyl.calc_mass(max(t_a + t_b, t_buckl))
        err = abs(M - M_i)
        M_i = M

    print(f'thickness stress: {(t_a + t_b):.4f} m\nthickness buckling: {t_buckl:.4f} m')
    print(f'MASS TOWER: {M / 1000:.4f} tonnes\n')



    # ---------------------------------------------------------------------------------
    print('\nPLATFORM')
    d = 150  # [m]
    r = d / 2
    h = 25 + 60


    # M = 1000e3
    M_i = 1000e3
    err = np.infty

    tower = Pylon(radius_outer=8, length=60 + 25, material=OS_Steel())
    M_platform = np.pi * tower.mat.rho * 0.5 * (r ** 2 - (r - 1) ** 2)

    while err > 1000:
        W = (M_truss + M_RNA + M_i + M_platform) * 9.81
        M_T = T_rated * (350 / 2)
        Fplx = T_rated / 4
        Fplz = (2 * M_T + r * W) / (4 * r)
        M_max = Fplx * h

        t_a = tower.axial_thickness(axial_force=Fplz)
        t_b = tower.bending_thickness(point_force=Fplx, position=h)
        t_buckl = tower.buckling(axial_load=Fplz)
        M = pyl.calc_mass(max(t_a + t_b, t_buckl))

        err = abs(M - M_i)
        M_i = M
        # print('iteration', err)

    print(f'thickness stress: {(t_a + t_b):.4f} m\nthickness buckling: {t_buckl:.4f} m')
    print(f'MASS TOWERS+PLATFORM: {(M/1000*4 + M_platform/1000):.4f} tonnes')
    print(f'MASS PLATFORM: {(M_platform / 1000):.4f} tonnes')
    print(f'MASS TOWERS: {(M/1000*4):.4f} tonnes')


    # ---------------------------------------------------------------------------------
    print('\nBRANCHING')
    L_beam = 3 * D_rot
    W_rna = M_RNA / 33

    a_arr = np.array([D_rot, 2 * D_rot, 3 * D_rot])
    Ts = np.ones(3) * T_perR
    Ws = W_rna * np.ones(3)
    thickness = 0.2

    branch = Pylon(radius_outer=1.5, length=L_beam, material=OS_Steel())
    delta_weight = branch.apply_loads(load_arr=Ws, pos_arr=a_arr, thickness=thickness)
    delta_thrust = branch.apply_loads(load_arr=Ts, pos_arr=a_arr, thickness=thickness)
    print(f'deflection thrust: {delta_thrust:.4f} m\ndeflection weight: {delta_weight:.4f} m')
    M = branch.calc_mass(thickness)
    print(M / 1000 * 12)

