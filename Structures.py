import numpy as np
import matplotlib.pyplot as plt

R = 29.6 # [m]
D = 2*R
N_rot  = 33


V_infty = 7
P_rated = 30e6
T_rated = P_rated/V_infty
T_perR = T_rated/N_rot


L_beam = 3*D
W_rna = 38790.46184*9.81


class OS_Steel():
    def __init__(self):
        self.E = 205e9
        self.rho = 7800
        self.sig_y = 470e6


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




class Pylon:
    def __init__(self, radius_outer, length, material):
        self.r_out = radius_outer
        self.L = length
        self.mat = material

    def calc_area(self):
        return 2*np.pi*self.r_out

    def calc_I(self):
        return np.pi*self.r_out**3

    def axial_thickness(self, axial_force):
        t_a = axial_force / (self.mat.sig_y * 2 * np.pi * self.r_out)
        return t_a

    def bending_thickness(self, point_force, position):
        M = position * point_force
        t_b = M / (np.pi * self.r_out ** 2 * self.mat.sig_y)
        return t_b

    def calc_mass(self):
        M = self.calc_area()*self.L*self.mat.rho
        return M

    def buckling(self, axial_load, n=0.25):
        t = self.L**2 * axial_load / (np.pi**3 * n * self.mat.E * self.r_out**3)
        return t

if __name__ == "__main__":
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