import numpy as np
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

'''ASSUMPTIONS:
=> thin-walled: D/t >= 10
=> slenderness ratio sufficiently high (col. buckling)
=> F_drag << F_thrust
=> F_lift << F_weight
=> D/t limits manufacturability
=> foundation mass not included 
=> safety factor accounts for limited load case [cite]
=> worst case scenario: thrust force is fully transferred as shear to top of tower
=> worst case: selft weight contributed to buckling
'''


class Steel():
    def __init__(self):
        'source: [cite]'
        self.E = 190e9      #[Pa]
        self.rho = 7800     #[kg/m^3]
        self.sig_y = 340e6  #[Pa]


class Sizing(ABC):
    def __init__(self, length, material):
        '''
        :param length: cylinder length [m]
        :param material: materials object
        '''
        self.L = length
        self.mat = material
        self.g = 9.80665

    @abstractmethod
    def axial_loading_function(self):
        pass

    @abstractmethod
    def bend_loading_function(self):
        pass

    @abstractmethod
    def calc_total_mass(self, M_tower):
        pass

    def calc_area(self, t, R_out):
        '''
        :param t: wall thickness [m]
        :return: CS area
        '''
        return 2 * np.pi * R_out * t

    def calc_I(self, t, R_out):
        '''
        :param t: thickness [m]
        :param R_out: outer radius [m]
        :return: moment of inertia (thin walled)
        '''
        return np.pi * R_out ** 3 * t

    def calc_mass(self, t, R_out):
        '''
        :param t: wall thickness [m]
        :param R_out: outer radius
        :return: mass [kg]
        '''
        M = self.calc_area(t=t, R_out=R_out) * self.L * self.mat.rho
        return M

    def column_buckling(self, F_internal, Length, rt, gamma_m_buckling = 1.2, n=0.25):
        '''
        :param Length: column length buckling analysis (i.e. load application point)
        :param n: end fixity
        :return: min. radius R for column buckling
        '''
        D = (8*abs(F_internal)*Length**2*gamma_m_buckling/(np.pi**3 *self.mat.E*rt*n))**(1/4)
        R=D/2
        return R

    @staticmethod
    def stress_radius(coefs, verbose: bool):
        '''
        :param coefs:
        :param verbose:
        :return:
        '''
        roots = np.roots(p=coefs)
        R_ult_stress = np.real(roots)[np.imag(roots) == 0][0]
        if verbose:
            print(f'Im(roots)={np.imag(roots)}')
        return R_ult_stress

    def sizing_analysis(self, rt=1/200, gamma_m = 1.1, gamma_m_buckling = 1.2, gamma_T=1.5, gamma_a=1.35, verbose=False):
        '''
        :param rt: t/D ratio
        :param gamma_m: material safety factor
        :param gamma_m_buckling: buckling safety factor
        :param gamma_T: thrust force safety factor
        :param gamma_a: axial safety factor
        :param verbose: print (bool)
        :return: sizing radius, sizing thickness
        '''
        if verbose:
            print(self)
        F_axial = self.axial_loading_function()*gamma_a
        Moment = self.bend_loading_function()*gamma_T
        coefs = np.array([self.mat.sig_y-self.L*self.mat.rho*self.g*gamma_m, 0, (-F_axial*gamma_m/(4*np.pi*rt)), (-Moment*gamma_m/(2*np.pi*rt))])

        R_ult_stress = float(self.stress_radius(coefs=coefs, verbose=verbose))
        M_tower = self.calc_mass(t=2*R_ult_stress*rt, R_out=R_ult_stress)
        R_col_buckling = self.column_buckling(F_internal=F_axial+self.g*M_tower, Length=self.L, rt=rt, gamma_m_buckling=gamma_m_buckling)

        R_final = max(R_ult_stress, R_col_buckling)

        M_tower = self.calc_mass(t=2*R_final*rt, R_out=R_final)
        M_final = self.calc_total_mass(M_tower=M_tower)
        #print(self.calc_area(t=2*R_final*rt, R_out=R_final))
        if verbose:
            print(f'D={2*R_final:.2f} R={R_final:.2f} R_str={R_ult_stress:.2f} R_buckl={R_col_buckling:.2f} '
                  f'M={M_final/1000:.2f} Mtower={M_tower/1000:.2f} t={rt*2*R_final:.4f}\n')
        return R_final, 2*R_final*rt, M_final


class SingleTower(Sizing):
    def __init__(self, length, material, M_truss, sum_M_RNA, F_T,):
        '''
        :param length: [m]
        :param material: object
        :param M_truss: [kg]
        :param sum_M_RNA: total RNA mass, all rotors [kg]
        :param F_T: average thrust force
        '''
        super().__init__(length=length, material=material)
        self.M_truss = M_truss
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T

    def __str__(self):
        return f'(i) Tower+Truss: L={self.L:.2f}'

    def axial_loading_function(self):
        '''
        :return: max internal axial force EXCLUDING SELF WEIGHT
        '''
        W_top = self.g * (self.M_truss + self.sum_M_RNA)
        return W_top

    def bend_loading_function(self):
        '''
        :return: max internal bending moment
        '''
        M = self.F_T * self.L
        return M

    def calc_total_mass(self, M_tower):
        '''
        :param M_tower: mass [kg] of a single tower
        :return: total concept mass
        '''
        return M_tower + self.M_truss


class Platform(Sizing):
    def __init__(self, length, material, M_truss, sum_M_RNA, F_T, r_platform, h_truss):
        '''
        :param length: [m]
        :param material: object
        :param M_truss: [kg]
        :param sum_M_RNA: total RNA mass, all rotors [kg]
        :param F_T: average total thrust force [N]
        '''
        super().__init__(length=length, material=material)
        self.M_truss = M_truss
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T
        self.r_platform = r_platform
        self.h_truss = h_truss
        self.M_platform = self.mass_platform()

    def __str__(self):
        return f'(ii) Platform+Truss: L={self.L:.2f}, r_plat={self.r_platform:.2f}, M_plat={self.M_platform/1000:.2f}'

    def mass_platform(self):
        tp = 0.5
        t_ring = 1
        return np.pi * self.mat.rho * tp * (self.r_platform ** 2 - (self.r_platform - t_ring) ** 2)

    def axial_loading_function(self):
        W = (self.M_truss + self.sum_M_RNA + self.M_platform) * 9.81
        M_T = self.F_T * self.h_truss / 2
        Fplz = (2 * M_T + self.r_platform * W) / (4 * self.r_platform)
        return Fplz

    def bend_loading_function(self):
        Fplx = self.F_T / 4
        M = Fplx*self.L
        return M

    def calc_total_mass(self, M_tower):
        #print(4*M_tower + self.M_truss + self.M_platform)
        return 4*M_tower + self.M_truss + self.M_platform


class Branches(Sizing):
    def __init__(self, length, material, sum_M_RNA, F_T, N_branches: int, N_rotors: int):
        '''
        :param length: [m]
        :param material: object
        :param M_truss: [kg]
        :param sum_M_RNA: total RNA mass, all rotors [kg]
        :param F_T: average thrust force
        '''
        super().__init__(length=length, material=material)
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T
        self.N_branches = N_branches
        self.N_rotors = N_rotors

    def axial_loading_function(self):
        '''
        :return: max internal axial force
        '''
        return 0

    def bend_loading_function(self):
        '''
        :return: max internal bending moment
        '''
        W_perR = self.sum_M_RNA/self.N_rotors
        T_perR = self.F_T/self.N_rotors
        M = 9/5*self.L*(W_perR+T_perR)
        return M

    def calc_total_mass(self, M_tower):
        '''
        :param M_tower: mass [kg] of a single tower
        :return: total concept mass
        '''
        return self.N_branches*M_tower


class BranchTower(Sizing):
    def __init__(self, length, material, sum_M_RNA, F_T, M_branches):
        '''
        :param length: [m]
        :param material: object
        :param M_truss: [kg]
        :param sum_M_RNA: total RNA mass, all rotors [kg]
        :param F_T: average thrust force
        '''
        super().__init__(length=length, material=material)
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T
        self.M_branches = M_branches

    def __str__(self):
        return f'(iii) Branching: L_tower={self.L:.2f}, M_branches={self.M_branches/1000:.2f}'

    def axial_loading_function(self):
        '''
        :return: max internal axial force
        '''
        return self.g*(self.M_branches + self.sum_M_RNA)

    def bend_loading_function(self):
        '''
        :return: max internal bending moment
        '''
        M = self.F_T*self.L
        return M

    def calc_total_mass(self, M_tower):
        '''
        :param M_tower: mass [kg] of a single tower
        :return: total concept mass
        '''
        return M_tower + self.M_branches


def statistics():
    references = {}


    'IEA 22 MW'
    dict = {}
    dict['hub_height'] = 170
    dict['tower_mass'] = 1574
    dict['monopile_mass'] = 2097
    dict['max_t'] = 0.091
    dict['min_t'] = 0.038
    dict['max_D'] = 10
    dict['min_D'] = 6
    dict['h_above_mud'] = 34+170
    dict['m_top'] = (821.2+120+82*3)*1000
    dict['P_rated'] = 22e6
    dict['v_rated'] = 11
    dict['T_rated'] = 2793e3

    references['IEA22MW'] = dict


    '5 MW reference'
    dict = {}
    dict['hub_height'] = 90
    dict['tower_mass'] = 347460/1000
    dict['monopile_mass'] = 0
    dict['max_t'] = .027
    dict['min_t'] = .019
    dict['max_D'] = 6
    dict['min_D'] = 3.87
    dict['h_above_mud'] = 87.6
    dict['P_rated'] = 18.45e6
    dict['v_rated'] = 11.4
    dict['T_rated'] = 800e3

    dict['m_top'] = 110e3 + 240e3
    references['5MWref'] = dict

    'IEA 15 MW'
    dict = {}
    dict['hub_height'] = 150
    dict['tower_mass'] = 860
    dict['monopile_mass'] = 1318
    dict['max_t'] = .055
    dict['min_t'] =.02
    dict['max_D'] = 10
    dict['min_D'] = 6.5
    dict['h_above_mud'] = 150+30
    dict['P_rated'] = 31.77#15e6
    dict['v_rated'] = 10.59
    dict['m_top'] = 1017e3
    dict['T_rated'] = 2.8e6
    references['IEA15MW'] = dict

    'DTU 8 MW'
    dict = {}
    dict['hub_height'] = 110
    dict['tower_mass'] = 558e3
    dict['monopile_mass'] = 0
    dict['max_t'] = .036
    dict['min_t'] = .022
    dict['max_D'] = 7.7
    dict['min_D'] = 5
    dict['h_above_mud'] = 110
    dict['P_rated'] = 34.2875e6
    dict['v_rated'] = 12.5
    dict['m_top'] = 90e3 + 285e3+35e3*3
    dict['T_rated'] = 2743e3
    references['DTU8MW'] = dict

    #print(references.keys())
    fig, axs = plt.subplots(1, 2, figsize=(10, 5),  layout='constrained')

    plus1 = []
    minus1 = []
    plus2 = []
    minus2 = []
    for key in references.keys():

        turbine = references[key]
        axs[0].plot(turbine['h_above_mud']*np.ones(2), [turbine['max_D'], turbine['min_D']], label=str(key))
        axs[1].plot(turbine['h_above_mud'] * np.ones(2), [turbine['max_t'], turbine['min_t']], label=str(key))

        single_tower = SingleTower(length=turbine['h_above_mud'], material=Steel(), M_truss=turbine['m_top'], sum_M_RNA=0, F_T=turbine['T_rated'], )
        R_single, t_single, _ = single_tower.sizing_analysis(verbose=False)
        D = 2 * R_single
        t = t_single
        axs[0].plot(turbine['h_above_mud'], D, marker='x')
        axs[1].plot(turbine['h_above_mud'], t, marker='x')
        plus1.append(abs(turbine['max_D'] - D)/D)
        minus1.append(abs(turbine['min_D'] - D)/D)

        plus2.append(abs(turbine['max_t'] - t) / t)
        minus2.append(abs(turbine['min_t'] - t) / t)

    #print(np.average(plus1)*100, np.average(minus1)*100)
    #print(np.average(plus2) * 100, np.average(minus2) * 100)
    axs[0].set_xlabel('cylinder height')
    axs[1].set_xlabel('cylinder height')

    axs[0].set_ylabel('Diameter [m]')
    axs[1].set_ylabel('thickness [m]')
    axs[0].legend()
    axs[1].legend()
    axs[0].set_title(f'uncertainty: +{np.average(plus1)*100:.2f}%, -{np.average(minus1)*100:.2f}%')
    axs[1].set_title(f'uncertainty: +{np.average(plus2) * 100:.2f}%, -{np.average(minus2) * 100:.2f}%')
    plt.show()

    uncertainty_range_D = [np.average(plus1), np.average(minus1)]
    uncertainty_range_t = [np.average(plus2) , np.average(minus2)]
    return uncertainty_range_D, uncertainty_range_t


if __name__ == '__main__':

    'UNCERTAINTY RANGES'
    uncertainty_range_D, uncertainty_range_t = statistics()

    'PARAMS:'
    R_rot = 29.6  # [m]
    D_rot = 2 * R_rot
    N_rot = 33

    V_infty = 10
    P_rated = 30e6
    T_rated = P_rated / V_infty
    T_perR = T_rated / N_rot

    M_truss = 2235959.595  # [kg]
    sum_RNA = 1280085.241  # [kg]

    max_depth = 60  # [m]
    frame_height = 350  # [m]
    clearance_height = 25  # [m]

    '(i) SINGLE TOWER'
    h_single = max_depth + clearance_height + frame_height / 2
    single_tower = SingleTower(length=h_single, material=Steel(), M_truss=M_truss, sum_M_RNA=sum_RNA, F_T=T_rated,)
    R_single, t_single, M_single = single_tower.sizing_analysis(verbose=True)
    print(M_single*(1+uncertainty_range_D[0])*(1+uncertainty_range_t[0])/1000, M_single*(1-uncertainty_range_D[1])*(1-uncertainty_range_t[1])/1000)

    '(ii) PLATFORM'
    h_platform = max_depth + clearance_height
    platform = Platform(length=h_platform, material=Steel(), M_truss=M_truss, sum_M_RNA=sum_RNA, F_T=T_rated, r_platform=50, h_truss=frame_height)
    _,_,M_plat_tot = platform.sizing_analysis(verbose=True)
    print(M_plat_tot * (1 + uncertainty_range_D[0]) * (1 + uncertainty_range_t[0]) / 1000,
          M_plat_tot * (1 - uncertainty_range_D[1]) * (1 - uncertainty_range_t[1]) / 1000)

    '(iii) BRANCHING'
    L_branch = 5 * R_rot
    h_branch_tow = max_depth + frame_height + clearance_height
    branches = Branches(length=L_branch, material=Steel(), sum_M_RNA=sum_RNA, F_T=T_rated, N_branches=12, N_rotors=N_rot)
    _, _, M_branches = branches.sizing_analysis()
    BrTower = BranchTower(length=h_branch_tow, material=Steel(), sum_M_RNA=sum_RNA, F_T=T_rated, M_branches=M_branches)
    _,_, M3 = BrTower.sizing_analysis(verbose=True)
    print(M3 * (1 + uncertainty_range_D[0]) * (1 + uncertainty_range_t[0]) / 1000,
          M3 * (1 - uncertainty_range_D[1]) * (1 - uncertainty_range_t[1]) / 1000)

    '''
    M_truss = 1000e3#2235959.595  # [kg]
    sum_RNA = 0*1280085.241  # [kg]
    P_rated = 15e6
    T_rated = P_rated / 10.59
    single_tower = Single_Tower(length=180, material=Steel(), M_truss=M_truss, sum_M_RNA=sum_RNA, F_T=T_rated, )
    single_tower.sizing_analysis(verbose=True)
    '''

    'DEFLECTION?'
    E = Steel().E
    print(f'\n(i) deflection, single = {T_rated * h_single ** 3 / (3 * E * np.pi * R_single ** 3 * t_single):.3f} [m]')
