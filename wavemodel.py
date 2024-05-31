import numpy as np
import matplotlib.pyplot as plt

class Wave():

    def __init__(self, lifetime, period, wavelength, water_depth, density, D, CM, CD, mu):
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
        self.mu = mu

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
    
    def compute_fluid_properties(self):

        Rn = self.wavelength/self.period*2*np.pi*self.D*self.rho/self.mu 
        KC = self.wavelength*2*np.pi/self.D
        beta = Rn/KC
        return beta, KC, Rn
    
    def compute_Morison(self, z, t):
        u = self.compute_water_velocity(z, t)
        u_dot = self.compute_water_acc(z, t)
        F = np.pi/4 * self.rho*self.CM*(self.D**2)*u_dot + 0.5*self.rho*self.D*self.CD*u*np.abs(u)
        return F



wave = Wave(lifetime=25, period = 5.2615, wavelength=121.1630, water_depth=20, density=1029, D = 8.5, CM = 1.7, CD= 0.6, mu = 1.8e-5)
z = np.arange(-wave.water_depth, 0, 0.1)
t = np.arange(0, 1000, 5)
print(wave.compute_fluid_properties())
Z, T = np.meshgrid(z,t)
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

F_eq = np.sum(0.1*F_dist_max)

Z_avg = np.sum(0.1*F_dist_max*z)/F_eq

M_base= np.sum 

print(F_eq, Z_avg+wave.water_depth, F_eq*(Z_avg+wave.water_depth))

'''fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter plot
ax.scatter(T, Z, F_dist)

# Set labels
ax.set_xlabel('Time [s]')
ax.set_ylabel('Depth [m]')
ax.set_zlabel('Force per unit length [N/m]')'''
plt.plot(F_dist_max, z)
plt.vlines(np.average(F_dist_max), ymin = -wave.water_depth, ymax=0)
plt.hlines(Z_avg, xmin = np.min(F_dist_max), xmax = np.max(F_dist_max))
plt.grid()
plt.xlabel('Force per unit length [N]')
plt.ylabel('Depth [m]')
# Show plot
plt.show()






