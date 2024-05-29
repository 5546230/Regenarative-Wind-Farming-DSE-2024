import numpy as np
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
from statistics_validation import statistics

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
    def __init__(self, length, material, M_truss, sum_M_RNA, F_T, F_L = 0, F_D = 0):
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
        self.F_L = F_L
        self.F_D = F_D

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

        M = (self.F_T) * self.L + self.F_L * 25  + self.F_D*(self.L+350/2)
        return M

    def calc_total_mass(self, M_tower):
        '''
        :param M_tower: mass [kg] of a single tower
        :return: total concept mass
        '''
        return M_tower + self.M_truss




if __name__ == '__main__':

    #'UNCERTAINTY RANGES'
    #uncertainty_range_D, uncertainty_range_t = statistics()

    'PARAMS:'
    R_rot = 29.6  # [m]
    D_rot = 2 * R_rot
    N_rot = 33

    V_infty = 10
    P_rated = 30e6
    T_rated = P_rated / V_infty
    T_perR = T_rated / N_rot

    M_truss = 7000e3#2235959.595  # [kg]
    sum_RNA = 1280085.241  # [kg]

    max_depth = 60  # [m]
    frame_height = 350  # [m]
    clearance_height = 25  # [m]

    '(i) SINGLE TOWER'
    h_single = 25+25+10#max_depth + clearance_height + frame_height / 2
    single_tower = SingleTower(length=h_single, material=Steel(), M_truss=M_truss, sum_M_RNA=sum_RNA, F_T=T_rated, F_L=3.4e6, F_D=300e3)
    R_single, t_single, M_single = single_tower.sizing_analysis(verbose=True)
    #print(M_single*(1+uncertainty_range_D[0])*(1+uncertainty_range_t[0])/1000, M_single*(1-uncertainty_range_D[1])*(1-uncertainty_range_t[1])/1000)


