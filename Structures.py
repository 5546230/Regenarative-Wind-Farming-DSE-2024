import numpy as np
import matplotlib.pyplot as plt

class Cyl_Geometry:
    def __init__(self, radius, length, thickness, material, thin: bool = False):
        self.r = radius
        self.l = length
        self.t = thickness
        self.mat = material
        self.thin = thin

    def calc_I(self):
        if not self.thin:
            return np.pi / 4 * self.r ** 4
        else:
            return np.pi / 4 * (self.r ** 4 - (self.r - self.t) ** 4)

    def calc_cs_area(self):
        A = np.pi * self.r**2
        if self.thin:
            A -= np.pi*(self.r-self.t)**2
        return A

    def calc_mass(self):
        M = self.mat.rho * self.calc_cs_area() * self.l
        return M



    def apply_loads(self, load_arr: np.array, pos_arr: np.array, thin: bool = False):
        EI = self.mat.E
        EI *= self.calc_I()

        delta_arr = load_arr * pos_arr * pos_arr / (6 * EI) * (3 * self.l - pos_arr)
        return sum(delta_arr)

    def calc_max_sig(self, P_b, P_a):
        M = -P_b*self.l
        max_sig_b = -M*self.r / self.calc_I()
        max_sig_a = P_a/self.calc_cs_area()

        return max_sig_b, max_sig_a


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
        return 2*np.pi*self.r_out*t

    def calc_I(self, t):
        '''
        :param t: thickness [m]
        :return: moment of inertia (thin walled)
        '''
        return np.pi*self.r_out**3*t

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
        M = self.calc_area(t=t)*self.L*self.mat.rho
        return M

    def buckling(self, axial_load, n=0.25):
        #print(self.mat.E)
        t = self.L**2 * axial_load / (np.pi**3 * n * self.mat.E * self.r_out**3)
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

    V_infty = 10
    P_rated = 30e6
    T_rated = P_rated / V_infty
    T_perR = T_rated / N_rot

    #W_rna = 38790.46184 * 9.81

    M_truss = 235959.595#*1.5
    M_RNA = 1280085.241

    print(f'MASS TRUSS: {M_truss/1000:.4f}\n')

    #---------------------------------------------------------------------------------
    print('TOWER')
    R_arr = np.linspace(5, 14, 20)
    mass_arr_1 = []
    for R in R_arr:
        #print(R)
        pyl = Pylon(radius_outer=R, length=60+25+350/2, material=OS_Steel())
        M_i = 1000e3
        err = np.infty
        while err>100:
            t_a = pyl.axial_thickness(axial_force=(M_truss+M_RNA+M_i)*9.81)
            t_b = pyl.bending_thickness(point_force=T_rated, position=60+25+.5*350)
            t_buckl = pyl.buckling(axial_load=(M_truss + M_RNA + 1*M_i)*9.81)
            M = pyl.calc_mass(max(t_a + t_b, t_buckl))
            err = abs(M-M_i)
            M_i = M


        #print(f'thickness stress: {(t_a + t_b):.4f} m\nthickness buckling: {t_buckl:.4f} m')
        if R / max(t_a + t_b, t_buckl) >= 10:
            mass_arr_1.append(M/1000 + M_truss/1000)
        else:
            mass_arr_1.append(0)


    plt.plot(R_arr, mass_arr_1, label='tower')
    #plt.show()
    #---------------------------------------------------------------------------------
    print('\nPLATFORM')
    d = 100 # [m]
    r = d/2
    h = 25+60

    R_arr = np.linspace(3, 14, 10)
    mass_arr_2 = []
    for R in R_arr:

        #M = 1000e3
        M_i = 1000e3
        err = np.infty

        tower = Pylon(radius_outer=R, length=60+25, material=OS_Steel())
        M_platform = np.pi * tower.mat.rho * 0.5 * (r ** 2 - (r - 1) ** 2)

        while err > 1000:
            W = (M_truss + M_RNA + (M_i*4) + M_platform) * 9.81
            M_T = T_rated * (350 / 2)
            Fplx = T_rated / 4
            Fplz = (2*M_T + r * W) / (4 * r)
            M_max = Fplx * h

            t_a = tower.axial_thickness(axial_force=Fplz)
            t_b = tower.bending_thickness(point_force=Fplx, position=h)
            t_buckl = tower.buckling(axial_load=Fplz)
            M = tower.calc_mass(max(t_a + t_b, t_buckl))

            err = abs(M - M_i)
            M_i = M
            #print('iteration', err)
        if R/max(t_a + t_b, t_buckl) >=10:
            mass_arr_2.append(M/1000*4 + M_platform/1000 + M_truss/1000)
        else:
            mass_arr_2.append(0)
        #print(f'thickness stress: {(t_a + t_b):.4f} m\nthickness buckling: {t_buckl:.4f} m')
        #print(f'MASS TOWERS+PLATFORM: {(M/1000*4 + M_platform/1000):.4f} tonnes')
        #print(f'MASS PLATFORM: {(M_platform / 1000):.4f} tonnes')
        #print(f'MASS TOWERS: {(M/1000*4):.4f} tonnes')

    plt.plot(R_arr, mass_arr_2, label=f'platform, R_platform = {r}')

    # ---------------------------------------------------------------------------------
    print('\nBRANCHING')
    L_beam = 3 * D_rot
    W_rna = M_RNA / 33

    a_arr = np.array([D_rot, 2 * D_rot, 3 * D_rot])
    Ts = np.ones(3)*T_perR
    Ws = W_rna * np.ones(3)
    thickness = 0.01

    M_max_thrust = 11*R_rot*T_perR

    R_arr = np.linspace(1, 10, 20)
    mass_arr_3 = []
    for R in R_arr:

        err = np.infty

        branch = Pylon(radius_outer=R, length=L_beam, material=OS_Steel())
        t_b = branch.bending_thickness(point_force=T_perR, position=11*R_rot)

        #t_buckl = tower.buckling(axial_load=Fplz)
        M = branch.calc_mass(t_b)
        if R / t_b >= 10:
            mass_arr_3.append(12*M/1000)

        else:
            mass_arr_3.append(0)


    plt.plot(R_arr, np.array(mass_arr_3), label=f'branching (only branches)')
    plt.legend()
    plt.xlabel('pylon Radius [m]')
    plt.ylabel('structural mass [tonnes]')
    plt.show()

    # ---------------------------------------------------------------------------------
    print('\nDTU TURBINE')
    L_beam = 3 * D_rot
    W_rna = M_RNA / 33
    h_hub = 150+30
    M_RNA = 1017e3

    T_rated = 15e6/V_infty

    R_arr = np.linspace(4, 10, 10)
    mass_arr_dtu = []
    for R in R_arr:

        # M = 1000e3
        M_i = 1000e3
        err = np.infty

        tower = Pylon(radius_outer=R, length=h_hub, material=OS_Steel())

        while err > 1000:

            t_a = tower.axial_thickness(axial_force=M_RNA*9.81)
            t_b = tower.bending_thickness(point_force=T_rated, position=h_hub)
            t_buckl = tower.buckling(axial_load=M_RNA)
            M = pyl.calc_mass(max(t_a + t_b, t_buckl))

            err = abs(M - M_i)
            M_i = M
            # print('iteration', err)
        if R / max(t_a + t_b, t_buckl) >= 10:
            mass_arr_dtu.append(M / 1000)
            #print(max(t_a + t_b, t_buckl))
        else:
            mass_arr_dtu.append(0)


    plt.plot(R_arr, np.array(mass_arr_dtu), label=f'DTU')
    plt.legend()
    plt.xlabel('member Radius [m]')
    plt.ylabel('structural mass [tonnes]')


    '''
    
    a_arr = np.array([D, 2 * D, 3 * D])
    Ts = np.array([T_perR, T_perR, T_perR])
    Ws = W_rna * np.ones(3)
    TEST = Cyl_Geometry(radius=2, length=L_beam, thickness=0.2, material=OS_Steel(), thin=True)
    print(f'Max tip deflection (thrust): {TEST.apply_loads(load_arr=Ts, pos_arr=a_arr)}')
    print(f'Max tip deflection (weight): {TEST.apply_loads(load_arr=Ws, pos_arr=a_arr)}')
    print(TEST.calc_mass() / 1000 * 12, 'tonnes')


    print()
    monopile_tower = Cyl_Geometry(radius=5, length = 300+60, thickness=0.3, material=OS_Steel(), thin=True)
    print(f'Max tip deflection (thrust): {monopile_tower.apply_loads(load_arr=np.array([T_rated]), pos_arr=np.array([60+175]))}')
    print(monopile_tower.calc_mass() / 1000,'tonnes' )


    print()
    single_tower = Cyl_Geometry(radius=5, length = 25 + 60, thickness=0.3, material=OS_Steel(), thin=True)
    print(f'Max tip deflection (thrust): {single_tower.apply_loads(load_arr=np.array([T_rated/3]), pos_arr=np.array([60 + 25]))}')
    print(single_tower.calc_mass() / 1000*3, 'tonnes')

    Pa = np.linspace(0, 3e9, 100)
    Pb = np.linspace(0, 1e9, 100)

    b, a = single_tower.calc_max_sig(P_b = Pb, P_a=Pa)

    plt.plot(Pa,a, label='axial')
    plt.plot(Pb,b, label='bending')
    plt.plot(Pa, np.ones_like(a)*single_tower.mat.sig_y)
    plt.legend()
    plt.show()
    '''