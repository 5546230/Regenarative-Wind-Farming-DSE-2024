import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


class Tower:
    def __init__(self, outer_radius, thickness, mass):
        self.R_o = outer_radius  
        self.t = thickness 
        self.m = mass  
        self.R_i = self.R_o - self.t 
        self.I_z = self.calculate_moment_of_inertia()

    def calculate_moment_of_inertia(self):
        
        I_z = 0.5*self.m*(self.R_i**2 + self.R_o**2)
        return I_z

class Rotor:
    def __init__(self, radius, n_blades, m_nacelle,):
        self.R = radius
        self.B = n_blades  
        self.m = m_nacelle + self.calculate_blades_mass()
    
    def calculate_blades_mass(self): 
        mass_blade = 2.6*1000
        total_mass_blades = self.B *mass_blade
        return total_mass_blades
    
class Array:
    def __init__(self, rotor_params, tower_params, array_params, V=6):
        self.g = 9.81  # gravity in m/s^2
        self.n_rows = array_params['n_rows']
        self.n_cols = array_params['n_cols']
        self.offset = array_params['offset']
        self.m_frame = array_params['m_frame']
        self.tower = Tower( outer_radius=tower_params['outer_radius'], 
                           thickness=tower_params['thickness'], mass=tower_params['mass'])
        self.rotor = Rotor(radius=rotor_params['radius'], n_blades=rotor_params['n_blades'],m_nacelle=rotor_params['m_nacelle'])
        
        self.spacing = 2 * self.rotor.R * (1+ self.offset)
        self.grid = self.create_grid()
        self.J =  self.calculate_I_z_array()
        self.q = 0.5 * 1.225 * V **2 

    def calculate_I_z_array(self):
        x_positions = self.grid[:, :, 0]
        m_per_turbine = self.rotor.m + (self.m_frame / (self.n_rows * 2* self.n_cols))
        I_z = np.sum(m_per_turbine * (x_positions**2))
        return I_z

    def create_grid(self):
        x0 = self.tower.R_o + self.offset + self.rotor.R
        z_coords = np.arange(self.n_rows) * self.spacing
        x_coords = np.concatenate([-x0 - np.arange(self.n_cols) * self.spacing, x0 + np.arange(self.n_cols) * self.spacing])

        grid = np.zeros((self.n_rows, 2 * self.n_cols, 3))
        for i in range(self.n_rows):
            grid[i, :, 0] = x_coords
            grid[i, :, 2] = z_coords[i]

        return grid
    
    def calc_thrust(self, c_T, phi):
        mag = c_T * self.q * (self.rotor.R**2 * np.pi)* np.cos(phi)**2
        F_t = np.zeros((self.n_rows, 2 * self.n_cols, 3))
        if phi > 0:
            # thrust on the negative x-side to create negative moment (yaw left)
            F_t[:, :self.n_cols, 1] = mag
            F_t[:, self.n_cols:, 1] = 0
        else:
            # thrust on the positive x-side to create positive moment (yaw right)
            F_t[:, self.n_cols:, 1] = mag
            F_t[:, :self.n_cols, 1] = 0
        return F_t

    def calc_moments(self,F_t):
        #y awing + tilting moment
        M = np.cross(self.grid.reshape(-1, 3), F_t.reshape(-1, 3)).reshape(self.n_rows, 2 * self.n_cols, 3)

        weight_per_unit = -self.g * (self.m_frame / (self.n_rows * 2 * self.n_cols) + self.rotor.m)
        W = np.zeros((self.n_rows, 2 * self.n_cols, 3))
        W[:, :, 2] = weight_per_unit  # weight acting downwards
        W = np.sum(W)
        
        Mx = np.sum(M[:, : , 0])
        Mz = np.sum(M[:, : , 2])
        #formula from kaydon catalog, idk if M_k is tilting of yawing moment
        
        return Mx, Mz, W
    
    def simulate_dynamics(self, c_T, total_time=300, dt=0.1):
        mu = 0.006
        phi = -0.17 # initial rad
        dphi = 0.000# initial rad/s
        times = np.arange(0, total_time, dt)
        phis = np.zeros_like(times)
        dphis = np.zeros_like(times)
        Mzs = np.zeros_like(times)
        Mws = np.zeros_like(times)

        yaw_motion_threshold = 0.000  # rad/s, adjust based on realistic system sensitivity

        for i, t in enumerate(times):
            T = self.calc_thrust(c_T, phi)
            Mx, Mz, W = self.calc_moments(T)
            ###motor test
            if i< 300:
                Mz = 2.9*15 * 1e6
            else:
                Mz=0
            ######
            G = 0  # Simplified assumption of no gyroscopic effects

            # Calculate dynamic friction moment only if there is significant motion
            if abs(dphi) > yaw_motion_threshold:
                Mw = -np.sign(dphi) * mu * (4.4 * (np.abs(Mx) + G) * 1.356 + 
                                        (np.abs(W) + self.tower.m * self.g) * 4.44822 * 2 * self.tower.R_i * 0.3048 +
                                        2.2 * np.abs(np.sum(T)) * 4.44822 * 2 * self.tower.R_i * 0.3048) / 2
            else:
                Mw = 0

            # eliminate numerical errors
            d2phi = (Mz + Mw) / self.J

            dphi += d2phi * dt
            phi += dphi * dt

            phis[i] = phi
            dphis[i] = dphi
            Mzs[i] = Mz
            Mws[i] = Mw
            if i == 0 or i % 500 == 0:
                print(Mz)

        return times, phis, dphis, Mzs, Mws


    def plot_Ms(self, times, Mzs, Mws):
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.plot(times, Mzs, label='M_z')
        plt.xlabel('Time (s)')
        plt.ylabel(r'$M_z$')
        plt.title('Moment over Time')
        plt.grid(True)
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.plot(times, Mws, label=r'$M_w$')
        plt.xlabel('Time (s)')
        plt.ylabel(r'$M_w$')
        plt.title('Friction over Time')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()

    def plot_results(self, times, phis, dphis):
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.plot(times, np.rad2deg(phis), label='Angle (phi)')
        plt.xlabel('Time (s)')
        plt.ylabel(r'$\phi$')
        plt.title('Angle over Time')
        plt.grid(True)
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.plot(times, np.rad2deg(dphis), label=r'$\frac{d\phi}{dt}$')
        plt.xlabel('Time (s)')
        plt.ylabel(r'$\frac{d\phi}{dt}$')
        plt.title('Angular Velocity over Time')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()

    def plot_array(self):
        x_coords = self.grid[:, :, 0].flatten()
        z_coords = self.grid[:, :, 2].flatten()
        
        plt.figure(figsize=(10, 5))
        plt.scatter(x_coords, z_coords, c='blue', marker='o', label='Turbines')
        plt.xlabel('X Coordinate (m)')
        plt.ylabel('Z Coordinate (m)')
        plt.title('Wind Turbine Array Layout')
        plt.axhline(0, color='grey', lw=0.5)  # Add y-axis line for reference
        plt.axvline(0, color='grey', lw=0.5)  # Add x-axis line for reference
        plt.grid(True)
        plt.legend()
        plt.show()
       




rotor_params = {
    'radius': 29.7,
    'n_blades': 3,
    'm_nacelle': 615*1000
}
array_params ={
    'n_rows': 4,
    'n_cols': 4,  # example density in kg/m^3
    'offset': 0.05,
    'm_frame': 2200*1000
}
tower_params = {
    'outer_radius': 8,
    'thickness': 0.3,
    'mass': 2200*1000  # Example mass in kg
}
array = Array(rotor_params,tower_params, array_params)
#array.plot_array()
c_T = 0.8
T = array.calc_thrust(c_T, phi=0.1)
times, phis, dphis, Mzs, Mws = array.simulate_dynamics(c_T)
array.plot_results(times, phis, dphis)
array.plot_Ms(times, Mzs, Mws)
print(array.J)
array.plot_array()
