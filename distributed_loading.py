import numpy as np
import matplotlib.pyplot as plt


class Wind_Model:
    def __init__(self, Vwr, zr, z0, alpha_shear,):
        self.Vwr = Vwr
        self.zr = zr
        self.z0 = z0
        self.alpha_shear = alpha_shear

    def log_shear_profile(self, zs, derivative: bool = False):
        '''
        :param z: height MSL [m]
        :param t: time [s]
        :return: velocity profile u(z, t), logarithmic power law
        '''
        Vw_z = self.Vwr * np.log(zs/self.z0) / np.log(self.zr/self.z0)
        if derivative:
            return self.Vwr  / np.log(self.zr/self.z0) / zs
        else:
            return Vw_z

    def pow_law_shear_profile(self, zs):
        '''
        :param z: height MSL [m]
        :param t: time [s]
        :return: velocity profile u(z, t), logarithmic power law
        '''
        Vw_z = self.Vwr * (zs/self.zr)**self.alpha_shear

        return Vw_z

    @staticmethod
    def plot_profile(z, V, label:str = None, show = True):
        plt.plot(V, z, label=label)
        plt.ylabel('height, z [m]')
        plt.xlabel('Velocity, V(z) [m/s]')
        if show:
            plt.legend()
            plt.show()


class Distr_Loading:
    def __init__(self, wind_model, wave_model):
        self.wave_model = wave_model
        self.wind_model = wind_model


    def morison_equation(self, z, f_CM: callable, f_CD: callable, rho, f_Dz:callable):

        f_uzt = self.wave_model.log_shear_profile(zs=z)
        f_uPrimezt = self.wave_model.log_shear_profile(zs=z, derivative=True)
        F_distr_inertial = (np.pi/4) * rho * f_CM(zs=z) * f_Dz(zs=z) * f_uPrimezt
        F_distr_drag = 0.5 * rho * f_CD(zs=z) * f_uzt * np.abs(f_uzt)

        f_Ft = F_distr_inertial + F_distr_drag
        return f_Ft



if __name__ == '__main__':
    wm = Wind_Model(Vwr=7, zr=10, z0=0.001, alpha_shear=0.12)
    z = np.linspace(1, 1000, 100)
    Vz_log = wm.log_shear_profile(zs=z,)
    Vz_log_prime = wm.log_shear_profile(zs=z, derivative=True)
    Vz_PL = wm.pow_law_shear_profile(zs=z,)
    wm.plot_profile(z=z, V=Vz_log, show=False, label='log')
    wm.plot_profile(z=z, V=Vz_PL, label='power law')
    wm.plot_profile(z=z, V=Vz_log_prime, show=True, label='log prime')

    loading_model = Distr_Loading(wind_model=wm, wave_model=wm)
    wave_distr = loading_model.morison_equation(z = z, rho=1.225, f_CM=lambda zs: 0.5, f_CD=lambda zs: 0.5, f_Dz=lambda zs:10)
    print(wave_distr.shape)