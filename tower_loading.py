import numpy as np
import matplotlib.pyplot as plt
from wavemodel import Wave
from drag_calc import Drag

class Loads():

    def __init__(self, tower_mass, truss_mass, nacelle_mass, truss_h, truss_depth, w_clearance, D_grid):

        
        self.tower_mass = tower_mass
        self.truss_mass = truss_mass
        self.n_mass = nacelle_mass
        self.T_rated = 3.945e6
        self.D_truss = D_grid
        self.L_afc = 5e6
        self.D_afc = 5e5
        self.h = truss_h/2+w_clearance
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

    def __init__(self, t, w_clearance, M_applied, F_applied, w_depth, mat_E, mat_yield):

        
        self.t = t
        self.rho  = 7850
        self.L = w_clearance+w_depth
        self.M = M_applied
        self.F = F_applied
        self.E = mat_E
        self.sy = mat_yield

    def calc_area(self, D):

        return np.pi*D*self.t
    
    def calc_mass(self, D):
        return self.calc_area(D)*0.1*self.rho
    
    def calc_Ixx(self, D):
        return np.pi*((D/2)**3)*self.t
    
    def calc_comp_stress_bending(self, D):
        return self.M*(D/2)/self.calc_Ixx(D)
    
    def calc_comp_stress_f(self, D):
        return self.F/self.calc_area(D)
    
    def calc_comp_stress(self, D):
        return self.calc_comp_stress_bending(D)+self.calc_comp_stress_f(D)
    
    def comp_buckling_stress(self,D):
        return 0.25*self.E*self.calc_Ixx(D)*(np.pi**2)/(0.1**2)/self.calc_area(D)
    
class Internal_loads():
    def __init__(self, F_wave, D_truss, F_comp, M_truss, M_thrust, M_afc, thrust, water_depth, z):
        self.F_wave= F_wave
        self.D_truss = D_truss
        self.D_afc = 5e5
        self.F_z = F_comp
        self.M_truss = M_truss
        self.M_th = M_thrust
        self.M_afc = M_afc
        self.T = thrust
        self.w_depth = water_depth
        self.z = z

    def int_axial(self):

        return self.F_z
    
    def int_shear(self):

        Vy = self.D_truss + self.T + np.cumsum(self.F_wave[::-1])+self.D_afc
        return Vy[::-1]
    
    def int_moment(self):

        M_x = self.M_afc+self.M_th+self.M_truss+(self.T+self.D_truss+self.D_afc)*(self.w_depth+z)+np.cumsum(self.F_wave[::-1]*0.1*(self.w_depth+z))
        return M_x[::-1]
    


    

if __name__ == "__main__":

    wave = Wave(lifetime=25, period = 5.2615, wavelength=121.1630, water_depth=20, density=1029, D = 8.5, CM = 1.7, CD= 0.6, mu = 1.3e-3)
    drag = Drag(V =13.5, rho=1.225, D_truss=1)
    D_grid, _ = drag.compute_CD_front()
    print(D_grid)
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

    
    load = Loads(tower_mass= 2235959.595, truss_mass=7000000, nacelle_mass=2602578.544, truss_depth=25, truss_h= 215, w_clearance=25, D_grid=D_grid)
    
    M_afc, M_truss, M_thrust = load.moments()

    F_comp = load.comp_force()
    
    print(M_base+M_afc+M_thrust+M_truss)
    print(F_comp*1e-6)
    
    int_loads = Internal_loads(F_wave=F_dist_max, D_truss = D_grid, F_comp = F_comp, M_truss = M_truss, M_thrust=M_thrust,M_afc=M_afc+67e6, thrust = 3.945e6, water_depth= 25, z = z)
    Vy = int_loads.int_shear()
    Mx = int_loads.int_moment()
    
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(Vy, z)
    ax1.set_title("Internal Shear Force [N]")
    ax1.set_xlabel(r'Shear Force $V_y$ [N]')
    ax1.set_ylabel(r'Position along the tower [m]')
    ax2.plot(Mx, z)
    ax2.set_title("Internal Bending Moment [Nm]")
    ax2.set_xlabel(r'Bending Moment $M_x$ [Nm]')
    ax2.set_ylabel(r'Position along the tower [m]')
    plt.show()  

    tower = Tower(t = 0.05, w_clearance=25, M_applied=1.4*Mx[0], F_applied=1.4*F_comp, w_depth=wave.water_depth, mat_E=190e9, mat_yield = 340e6)
    D_dist  = np.ones(len(z))*8
    
    t_it = 0.05
    mass = np.sum(tower.calc_mass(D_dist))
    print(f'Original mass: {mass}')

    for i in range(len(D_dist)):
        D_it1 = D_dist[i]
        tower_it = Tower(t = 0.05, w_clearance=25, M_applied=1.4*Mx[i], F_applied=1.4*F_comp, w_depth=wave.water_depth, mat_E=190e9, mat_yield = 340e6)
        tower_stress = tower_it.calc_comp_stress(D=D_dist[i])
        buckling = tower_it.comp_buckling_stress(D=D_dist[i])
        mass = np.sum(tower_it.calc_mass(D_dist[i]))

        print(i/len(D_dist))
        while tower_stress > np.minimum(buckling, tower_it.sy):
            
            D_it1 += 0.01
            tower_stress = tower_it.calc_comp_stress(D_it1)
            
            #buckling = tower_it.comp_buckling_stress(D_it1)
            
            '''wave1 = Wave(lifetime=25, period = 5.2615, wavelength=121.1630, water_depth=20, density=1029, D = D_it1, CM = 1.7, CD= 0.6, mu = 1.3e-3)

            z = np.arange(-wave1.water_depth, 0, 0.1)
            t = np.arange(0, 1000, 5)
            F_dist=[]
            avg=0
            avg_lst=[]
            for j in t:

                F = wave1.compute_Morison(z, j)
                F_dist.append(F)
                avg = np.average(F)
                avg_lst.append(avg)

            F_dist = np.array(F_dist)
            avg = np.array(avg_lst)
            max_idx = np.argmax(avg)

            F_dist_max = F_dist[max_idx]

            
            #t_it += 0.001'''
           
            


            
            
        D_dist[i] = D_it1

    coef = np.polyfit(D_dist/2, z, deg = 3)
    z_fit1 = np.poly1d(coef)
    z_fitted1 = z_fit1(D_dist/2)

    coef2 = np.polyfit(-D_dist/2, z, deg = 3)
    z_fit2 = np.poly1d(coef2)
    z_fitted2 = z_fit2(-D_dist/2)
        
    #plt.plot(-D_dist/2, z, color = 'b')
    #plt.plot(D_dist/2, z, color = 'b')
    plt.plot(D_dist/2, z_fitted1, color = 'g')
    plt.plot(-D_dist/2, z_fitted2, color = 'g')
    plt.title("Monopile Sizing, thickness 5mm")
    plt.xlabel(r"Monopile Diameter [m] ")
    plt.ylabel(r'Position along the monopile [m]')


    
    plt.axis('equal')
    plt.show()

    mass = np.sum(tower.calc_mass(D_dist))
    mass_inop = tower.calc_mass(D_dist[0])*len(D_dist)
    print(f'Optimised mass: {mass}')
    print(f'Unoptimised mass: {mass_inop}')
    print(f'Mass savings: {mass_inop-mass}')
    print(D_dist[0])
    print(D_dist[-1])

            

            
                
            
    

    



