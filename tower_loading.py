import numpy as np
import matplotlib.pyplot as plt
from wavemodel import Wave
from drag_calc import D_element

class Loads():

    def __init__(self, tower_mass, truss_mass, nacelle_mass, truss_h, truss_depth, w_clearance, w_depth, M_base):

        self.M_wave = M_base
        self.tower_mass = tower_mass
        self.truss_mass = truss_mass
        self.n_mass = nacelle_mass
        self.T_rated = 3.975734e6
        self.D_truss = D_element*33
        self.L_afc = 5e6
        self.D_afc = 5e5
        self.h = truss_h/2+w_clearance+w_depth
        self.d = truss_depth
        self.g = 9.80665

    def moments(self):
        M_afc = self.L_afc*self.d + self.D_afc*self.h
        M_truss = self.D_truss * self.h
        M_thrust = self.T_rated * self.h

        return M_afc, M_truss, M_thrust
    def comp_force(self):

        return self.L_afc+(self.truss_mass+self.tower_mass+self.n_mass)*self.g
    
class Tower():

    def __init__(self, D, t, w_clearance, M_applied, F_applied, w_depth, mat_E, mat_yield):

        self.D = D
        self.t = t
        self.rho  = 7850
        self.L = w_clearance+w_depth
        self.M = M_applied
        self.F = F_applied
        self.E = mat_E
        self.sy = mat_yield

    def calc_area(self):

        return np.pi*self.D*self.t
    
    def calc_mass(self):
        return self.calc_area()*self.L*self.rho
    
    def calc_Ixx(self):
        return np.pi*((self.D/2)**3)*self.t
    
    def calc_comp_stress_bending(self):
        return self.M*(self.D/2)/self.calc_Ixx()
    
    def calc_comp_stress_f(self):
        return self.F/self.calc_area()
    
    def calc_comp_stress(self):
        return self.calc_comp_stress_bending()+self.calc_comp_stress_f()
    
    def comp_buckling_stress(self):
        return 0.25*self.E*self.calc_Ixx()*(np.pi**2)/(self.L**2)/self.calc_area()

    

if __name__ == "__main__":

    wave = Wave(lifetime=25, period = 5.2615, wavelength=121.1630, water_depth=20, density=1029, D = 8.5, CM = 1.7, CD= 0.6, mu = 1.3e-3)

    z = np.arange(-wave.water_depth, 0, 0.1)
    t = np.arange(0, 1000, 5)
    F_dist=[]
    avg=0
    avg_lst=[]
    for i in t:

        F = wave.compute_Morison(z, i)
        F_dist.append(F)
        avg = np.average(F)
        avg_lst.append(avg)

    F_dist = np.array(F_dist)
    avg = np.array(avg_lst)
    max_idx = np.argmax(avg)

    F_dist_max = F_dist[max_idx]

    M_base= np.sum((wave.water_depth+z)*F_dist_max*0.2)

    
    load = Loads(tower_mass= 2235959.595, truss_mass=5000000, nacelle_mass=2602578.544, truss_depth=25, truss_h= 215, w_clearance=25, w_depth=wave.water_depth, M_base=M_base)
    
    M_afc, M_truss, M_thrust = load.moments()

    F_comp = load.comp_force()

    print(M_base+M_afc+M_thrust+M_truss)
    print(F_comp*1e-6)
    tower = Tower(D=8.5, t = 0.05, w_clearance=25, M_applied=M_base+M_afc+M_thrust+M_truss, F_applied=F_comp, w_depth=wave.water_depth, mat_E=190e9, mat_yield = 340e6)
    
    tower_stress = tower.calc_comp_stress()
    buckling = tower.comp_buckling_stress()
    D_it = 8.5
    t_it = 0.05
    mass = tower.calc_mass()
    print(f'Original mass: {mass}')
    while tower_stress > np.minimum(buckling, tower.sy):
        D_it1 = D_it+0.001

        wave1 = Wave(lifetime=25, period = 5.2615, wavelength=121.1630, water_depth=20, density=1029, D = D_it1, CM = 1.7, CD= 0.6, mu = 1.3e-3)

        z = np.arange(-wave1.water_depth, 0, 0.1)
        t = np.arange(0, 1000, 5)
        F_dist=[]
        avg=0
        avg_lst=[]
        for i in t:

            F = wave1.compute_Morison(z, i)
            F_dist.append(F)
            avg = np.average(F)
            avg_lst.append(avg)

        F_dist = np.array(F_dist)
        avg = np.array(avg_lst)
        max_idx = np.argmax(avg)

        F_dist_max = F_dist[max_idx]

        M_base= np.sum((wave.water_depth+z)*F_dist_max*0.2)
        #t_it += 0.001
        load1 = Loads(tower_mass= mass, truss_mass=5000000, nacelle_mass=2602578.544, truss_depth=25, truss_h= 215, w_clearance=25, w_depth=wave.water_depth, M_base=M_base)
    
        M_afc, M_truss, M_thrust = load1.moments()

        F_comp = load1.comp_force()

        tower_it = Tower(D=D_it1, t = t_it, w_clearance=25, M_applied=M_base+M_afc+M_thrust+M_truss, F_applied=F_comp, w_depth=wave.water_depth, mat_E=190e9, mat_yield = 340e6)
        tower_stress = tower_it.calc_comp_stress()
        buckling = tower_it.comp_buckling_stress()
        mass = tower_it.calc_mass()
        D_it = D_it1

        

        
            
        
    print(f'Final mass:{mass}')
    print(f'Final Diameter: {D_it}')
    print(f'thickness: {t_it}')


    



