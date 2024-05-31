import numpy as np
import matplotlib.pyplot as plt

class Wave():

    def __init__(self, lifetime, period, wavelength, water_depth, density, D, CM, CD):
        '''params:
        lifetime in years'''
        self.lifetime = lifetime
        self.period = period
        self.wavelength = wavelength
        self.water_depth = water_depth
        self.rho = density
        self.D = D
        self.CM = CM
        self.CD = CD

    def compute_wave_properties(self):
        Hs = 0.479 * np.log(self.lifetime)+6.0626
        omega = 2*np.pi/self.period
        k = 2*np.pi/self.wavelength

        return Hs, omega, k



    def compute_water_velocity(self, z, t):
        Hs, omega, k = self.compute_wave_properties()
        return omega*Hs/2 * (np.cosh(k*(z+self.water_depth))/np.sinh(k*self.water_depth)) * np.cos(omega*t)
    
    def compute_water_acc(self, z, t):
        Hs, omega, k = self.compute_wave_properties()
        return -(omega**2)*Hs/2 * (np.cosh(k*(z+self.water_depth))/np.sinh(k*self.water_depth)) * np.sin(omega*t)
    
    def compute_Morison(self, z, t):
        u = self.compute_water_velocity(z, t)
        u_dot = self.compute_water_acc(z, t)
        F = np.pi/4 * self.rho*self.CM*(self.D**2)*u_dot + 0.5*self.rho*self.D*self.CD*u*np.abs(u)
        return F



wave = Wave(lifetime=25, period = 10, wavelength=30, water_depth=20, density=1029, D = 10, CM = 2, CD= 1.2)
z = np.arange(-wave.water_depth, 0, 0.2)
t = np.arange(0, 1000, 10)
F = wave.compute_Morison(z, t)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter plot
ax.scatter(F, z, t)

# Set labels
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Coordinate')

# Show plot
plt.show()




