'''
Improved way of power curve
'''

import numpy as np
from matplotlib import pyplot as plt

# +----------------------------+
# |          Inputs            |
# +----------------------------+
equivalent_radius = 170 # [m]
rho = 1.225 # [kg/m^3]
n_rotors = 33
P_RATED = 30e6 # [W]
V_RATED = 10 # [m/s]
TSR = 8 # [-]
CP = 0.46
cut_in = 3 # [m/s]
cut_off = 25 # [m/s]
#////////////////////////////////

# Calculate radius of multi rotor
radius = equivalent_radius / (n_rotors ** 0.5)

# Initialize Array
U_array = np.linspace(0, 27, 350)

# Calculate Power and Torque
P_array = 0.5 * rho * CP * U_array ** 3 * np.pi * radius ** 2 * n_rotors
Q_array = 0.5 * rho * CP * np.pi * radius ** 5 /(TSR ** 3) * (U_array * TSR / radius) ** 2

# Apply cut-in, rated and cut-off constraints
P_array[U_array < cut_in] = 0
Q_array[U_array < cut_in] = 0
Q_array[P_array > P_RATED] = Q_array[P_array > P_RATED][0]
P_array[P_array > P_RATED] = P_RATED
P_array[U_array > cut_off] = 0
Q_array[U_array > cut_off] = 0

P_array /= 1e6
Q_array /= 1e3

fig, axs = plt.subplots(1, 2)
fig.set_figheight(3)
fig.set_figwidth(10)

axs[0].plot(U_array, P_array)
axs[0].set_xlabel('Wind Speed [m/s]')
axs[0].set_xlim(0, 27)
axs[0].set_ylabel('Power Generated [MW]')
axs[0].grid(True)

axs[1].plot(U_array, Q_array)
axs[1].set_xlabel('Wind Speed [m/s]')
axs[1].set_xlim(0, 27)
axs[1].set_ylabel('Torque per Rotor [kNm]')
axs[1].grid(True)

plt.tight_layout()
plt.savefig('power_torque_curves.svg', format='svg')
plt.show()
