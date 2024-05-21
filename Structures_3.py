import numpy as np
import matplotlib.pyplot as plt

class Steel():
    def __init__(self):
        self.E = 190e9
        self.rho = 7800
        self.sig_y = 340e6

class Sizing:
    def __init__(self, length,  material):
        '''
        :param length: cylinder length [m]
        :param material: materials object
        '''
        self.L = length
        self.mat = material

    def calc_area(self, t, R_out):
        '''
        :param t: wall thickness [m]
        :return: CS area
        '''
        return 2*np.pi*R_out*t

    def calc_I(self, t, R_out):
        '''
        :param t: thickness [m]
        :param R_out: outer radius [m]
        :return: moment of inertia (thin walled)
        '''
        return np.pi*R_out**3*t

    def axial_thickness(self, axial_force, R_out):
        '''
        :param axial_force: applied axial force (N)
        :param R_out: outer radius [m]
        :return: wall thickness to sustain axial load within sigma_yield
        '''
        t_a = axial_force / (self.mat.sig_y * 2 * np.pi * R_out)
        return t_a

    def bending_thickness(self, point_force, position, R_out):
        '''
        :param point_force: point force magnitude [N]
        :param position: point of application along beam [m]
        :param R_out: outer radius [m]
        :return: t to withstand max bending moment at sigma_yield
        '''
        M = position * point_force
        t_b = M / (np.pi * R_out ** 2 * self.mat.sig_y)
        return t_b

    def calc_mass(self, t, R_out):
        '''
        :param t: wall thickness [m]
        :param R_out: outer radius
        :return: mass [kg]
        '''
        M = self.calc_area(t=t, R_out=R_out)*self.L*self.mat.rho
        return M

    def column_buckling(self, axial_load, R_out, n=0.25):
        '''
        :param axial_load: [N]
        :param R_out: [m]
        :param n: [-]
        :return: thickness
        '''
        t = self.L**2 * axial_load / (np.pi**3 * n * self.mat.E * R_out**3)
        return t


class Single_Tower(Sizing):
    def __init__(self, length, material, M_truss, sum_M_RNA, F_T, mass_0):
        '''
        :param length: [m]
        :param material: object
        :param M_truss: [kg]
        :param sum_M_RNA: total RNA mass, all rotors [kg]
        :param F_T: average thrust force
        :param mass_0: initial condition iterator
        '''
        super().__init__(length=length,  material=material)
        self.M_truss = M_truss
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T
        self.M_0 = mass_0

    def sizing_analysis(self, R_array: np.array, tol, SF, verbose = False):
        '''
        :param R_array: array of outer radius values
        :param tol: iter tolerance
        :param SF: safety factor
        :return: [-]
        '''
        mass_arr = []
        Dt_arr = []

        M_0 = self.M_0
        M_1 = M_0
        t_axial = t_bending = t_col_buckling = np.infty

        for R in R_array:
            err = np.infty

            while err>tol:
                weight_axial = (self.M_truss + self.sum_M_RNA + M_0) * 9.81
                t_axial = self.axial_thickness(axial_force=weight_axial*SF, R_out=R)
                t_bending = self.bending_thickness(point_force=self.F_T*SF, position=self.L, R_out=R)
                t_col_buckling = self.column_buckling(axial_load=(self.M_truss + self.sum_M_RNA +M_0*0)*9.81*SF, R_out=R)

                M_1 = self.calc_mass(t=max(t_axial + t_bending, t_col_buckling), R_out=R)
                err = abs(M_1 - M_0)
                M_0 = M_1

            Dt = 2*R/max(t_axial + t_bending, t_col_buckling)
            Dt_arr.append(Dt)
            mass_arr.append(M_1 + self.M_truss)

            if verbose:
                print(f'R = {R:.4f}, M = {(M_1 + self.M_truss) / 1000:.4f} t, t = {max(t_axial + t_bending, t_col_buckling):.4f}, D/t = {Dt:.4f}')
        return np.array(R_array), np.array(mass_arr), np.array(Dt_arr)


class Platform(Sizing):
    def __init__(self, length, material, M_truss, sum_M_RNA, F_T, mass_0, r_platform, h_truss):
        '''
        :param length: [m]
        :param material: obj
        :param M_truss: [kg]
        :param sum_M_RNA: [kg]
        :param F_T: [N]
        :param mass_0: [kg]
        :param r_platform: radius of platform [m]
        '''
        super().__init__(length=length,  material=material)
        self.M_truss = M_truss
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T
        self.M_0 = mass_0
        self.r_platform = r_platform
        self.M_platform = self.mass_platform()
        self.h_truss = h_truss

    def mass_platform(self):
        tp = 0.5
        t_ring = 1
        return  np.pi * self.mat.rho * tp * (self.r_platform ** 2 - (self.r_platform - t_ring) ** 2)

    def sizing_analysis(self, R_array: np.array, tol, SF, verbose=False):
        mass_arr = []
        Dt_arr = []

        M_0 = self.M_0
        M_1 = M_0
        t_axial = t_bending = t_col_buckling = np.infty
        'moment due to thrust force'
        M_T = self.F_T *self.h_truss/2

        for R in R_array:
            err = np.infty

            while err>tol:
                'weight (excluding tower self weight)'
                W = (self.M_truss + self.sum_M_RNA + self.M_platform) * 9.81

                'shear force'
                Fplx = T_rated / 4
                'sizing axial force'
                Fplz = (2 * M_T + self.r_platform * W) / (4 * self.r_platform)

                t_axial = self.axial_thickness(axial_force=Fplz*SF, R_out=R)
                t_bending = self.bending_thickness(point_force=Fplx*SF, position=self.L, R_out=R)
                t_col_buckling = self.column_buckling(axial_load=(Fplz +M_0*9.81*0)*SF, R_out=R)

                M_1 = self.calc_mass(t=max(t_axial + t_bending, t_col_buckling), R_out=R)
                err = abs(M_1 - M_0)
                M_0 = M_1

            Dt = 2*R/max(t_axial + t_bending, t_col_buckling)
            Dt_arr.append(Dt)
            #print(M_1*4/1000,  self.M_truss/1000,  self.M_platform/1000, (M_1*4 + self.M_truss + self.M_platform) / 1000)
            mass_arr.append(M_1*4 + self.M_truss + self.M_platform)

            if verbose:
                print(f'R = {R:.4f}, M = {(M_1*4 + self.M_truss + self.M_platform) / 1000:.4f} t, t = {max(t_axial + t_bending, t_col_buckling):.4f}, D/t = {Dt:.4f}')
        return np.array(R_array), np.array(mass_arr), np.array(Dt_arr)


class Branches(Sizing):
    def __init__(self, length, material, sum_M_RNA, F_T, N_branch, r_rotor, N_rotor, mass_0=1000):
        super().__init__(length=length,  material=material)
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T
        self.M_0 = mass_0
        self.N_branch = N_branch
        self.N_rotor = N_rotor
        self.r_rotor = r_rotor

        self.pos_array = np.array([r_rotor, 3*r_rotor, 5*r_rotor])
        self.T_array = F_T/N_rotor * np.ones(3)
        self.W_array = sum_M_RNA / N_rotor * np.ones(3)

    def _branch_sizing(self, R_branch: np.array, tol, SF,):
        M_0 = self.M_0
        M_1 = M_0

        P_max_weight = self.W_array[0]
        P_max_T = self.T_array[0]
        d =  9 * self.r_rotor

        t_bending_T = t_bending_weight = np.infty
        err = np.infty
        while err>tol:
            t_bending_weight = self.bending_thickness(point_force=P_max_weight*SF, position=d, R_out=R_branch)
            t_bending_T = self.bending_thickness(point_force=P_max_T*SF, position=d, R_out=R_branch)

            M_1 = self.calc_mass(t=max(t_bending_T, t_bending_weight), R_out=R_branch)
            err = abs(M_1 - M_0)
            M_0 = M_1

        t = max(t_bending_T, t_bending_weight)
        Dt = 2*R_branch/t

        deflection1 = self.deflection(P_arr=self.W_array, a_arr=self.pos_array, t=t, R_out = R_branch)
        deflection2 = self.deflection(P_arr=self.T_array, a_arr=self.pos_array, t=t, R_out=R_branch)
        deflection = max(deflection1, deflection2)

        return self.N_branch*M_1, Dt, max(t_bending_T, t_bending_weight), deflection

    def deflection(self, P_arr, a_arr, t, R_out):
        return np.sum(P_arr*a_arr**2 / (6*self.mat.E*self.calc_I(t=t, R_out=R_out))* (3*self.L-a_arr))


class Branch_Tower(Sizing):
    def __init__(self, length, material, sum_M_RNA, F_T, Branches: Branches, h_truss, mass_0=1000):
        super().__init__(length=length,  material=material)
        self.sum_M_RNA = sum_M_RNA
        self.F_T = F_T
        self.M_0 = mass_0
        self.Branches = Branches
        self.h_truss = h_truss

    def sizing_analysis(self, R_tow_array: np.array, R_branch_array, tol, SF, verbose=False):
        mass_arr = []
        Dt_arr_tow = []
        Dt_arr_branch = []

        M_0 = self.M_0
        M_1 = M_0
        t_axial = t_bending = t_col_buckling = branches_mass = np.infty
        Dt_branch=0
        for Rt in R_tow_array:
            for Rb in R_branch_array:
                err = np.infty

                while err>tol:
                    branches_mass, Dt_branch, tb, deflection = self.Branches._branch_sizing(R_branch = Rb, tol=100, SF=SF)
                    weight = ( self.sum_M_RNA + branches_mass) * 9.81

                    t_axial = self.axial_thickness(axial_force=weight * SF, R_out=Rt)
                    t_bending = self.bending_thickness(point_force=self.F_T * SF, position=self.L-self.h_truss/2, R_out=Rt)
                    t_col_buckling = self.column_buckling(axial_load=( self.sum_M_RNA + branches_mass) * SF * 9.81, R_out=Rt)

                    M_1 = self.calc_mass(t=max(t_axial + t_bending, t_col_buckling), R_out=Rt)
                    err = abs(M_1 - M_0)
                    M_0 = M_1

                t = max(t_axial + t_bending, t_col_buckling)
                Dt = 2*Rt/t
                Dt_arr_tow.append(Dt)
                Dt_arr_branch.append(Dt_branch)
                mass_arr.append(M_1 + branches_mass)

                if verbose and 305>Dt>295 and 305>Dt_branch>295:
                    print(f'Rt = {Rt:.4f}, Rb={Rb:.4f}, M = {(M_1 + branches_mass) / 1000:.4f} t, tt = {max(t_axial + t_bending, t_col_buckling):.4f}, '
                          f'tb={tb:.4f}, D/t_t = {Dt:.4f}, D/t_b = {Dt_branch:.4f}, delta={deflection:.4f}')
        return np.array(R_tow_array), np.array(R_branch_array), np.array(mass_arr), np.array(Dt_arr_tow), np.array(Dt_arr_branch)


def plot_analysis(Rs, ts, Ms, label: str, limit=300):
    fig, ax1 = plt.subplots()
    ax1.plot(Rs, Ms/1000, label=label)
    ax1.set_xlabel('radius (R) [m]')
    ax1.set_ylabel('mass  [tonnes]')

    ax2 = ax1.twinx()
    ax2.plot(Rs, ts, color='green', label='D/t')
    ax2.plot(Rs, np.ones_like(Rs)*limit, linestyle='--', color='green',label='D/t limit')
    ax2.set_ylabel('D/t')
    fig.legend()
    plt.show()



if __name__ == '__main__':
    R_rot = 29.6  # [m]
    D_rot = 2 * R_rot
    N_rot = 33

    V_infty = 10
    P_rated = 30e6
    T_rated = P_rated / V_infty
    T_perR = T_rated / N_rot

    M_truss = 2235959.595 #[kg]
    sum_RNA = 1280085.241 # [kg]
    max_depth = 60 # [m]
    frame_height = 350 # [m]
    clearance_height = 25 # [m]



    h_single = max_depth + clearance_height + frame_height/2
    single_tower = Single_Tower(length=h_single, material=Steel(), M_truss = M_truss, sum_M_RNA=sum_RNA, F_T=T_rated, mass_0=2000e3,)
    Rs1, Ms1, ts1 = single_tower.sizing_analysis(R_array=np.linspace(3., 10, 10), tol=50, SF=1, verbose=False)

    h_platform = max_depth + clearance_height
    platform = Platform(length=h_platform, material=Steel(), M_truss=M_truss, sum_M_RNA=sum_RNA, F_T=T_rated, mass_0=2000e3, r_platform=100, h_truss=frame_height)
    Rs2, Ms2, ts2 = platform.sizing_analysis(R_array=np.linspace(1., 3., 100), tol=50, SF=1, verbose=True)

    L_branch  = 5*R_rot
    h_branch_tow = max_depth + frame_height + clearance_height
    branches = Branches(length=L_branch, material=Steel(), sum_M_RNA=sum_RNA, F_T=T_rated, N_branch=12, r_rotor=R_rot, N_rotor=N_rot)
    branch_tower = Branch_Tower(length=h_branch_tow, material=Steel(), sum_M_RNA=sum_RNA, F_T = T_rated, Branches=branches, h_truss=frame_height)
    R_tower, R_branch, Ms3, ts_tower, ts_branch = branch_tower.sizing_analysis(R_branch_array=np.linspace(.5, 2., 100),
                                                                               R_tow_array=np.linspace(3., 8, 100),
                                                                               tol=50, SF=1, verbose=False)




    fig, ax1 = plt.subplots()
    CS = ax1.contourf(R_tower, R_branch, Ms3.reshape(R_tower.size, R_branch.size)/1000, levels=20)
    ax1.set_xlabel(r'$R_{tower}$')
    ax1.set_ylabel(r'$R_{branch}$')
    fig.colorbar(CS, label='Mass [tonnes]')
    ax1.contour(R_tower, R_branch, ts_tower.reshape(R_tower.size, R_branch.size), levels=np.linspace(295, 305, 1), colors='white')
    ax1.contour(R_tower, R_branch, ts_branch.reshape(R_tower.size, R_branch.size),  levels=np.linspace(295, 305, 1), colors='white')
    plt.show()

    plot_analysis(Rs = Rs1, ts=ts1, Ms=Ms1, label='tower+truss')
    plot_analysis(Rs=Rs2, ts=ts2, Ms=Ms2, label='platform+truss')



    '''
    fig, ax1 = plt.subplots()
    ax1.plot(Rs1, Ms1/1000, label='Tower+Truss')
    #ax1.plot(Rs2, Ms2/ 1000, label='Platform+Truss')
    #plt.plot(Rs3, Ms3 / 1000, label='Branches')
    ax1.set_xlabel('radius (R) [m]')
    ax1.set_ylabel('mass  [tonnes]')

    ax2 = ax1.twinx()
    ax2.plot(Rs1, ts1, color='green', label='D/t')
    ax2.plot(Rs1, np.ones_like(Rs1)*300, linestyle='--', color='green',label='D/t limit')
    fig.legend()
    plt.show()


    plt.plot(Rs1, ts1 , label='Tower+Truss')
    plt.plot(Rs2, ts2, label='Platform+Truss')
    plt.plot(Rs2, np.ones_like(Rs1)*300, linestyle='--', label='D/t limit')
    plt.xlabel('radius (R) [m]')
    plt.ylabel('D/t  [-]')
    plt.ylim([0, 500])
    plt.legend()
    plt.show()
    '''