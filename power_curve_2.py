'''
Improved way of power curve
'''

import numpy as np
from matplotlib import pyplot as plt

# +----------------------------+
# |          Inputs            |
# +----------------------------+
#equivalent_radius = 170 # [m]
rho = 1.225 # [kg/m^3]
n_rotors = 33
P_RATED = 30e6 # [W]
V_RATED = 10 # [m/s]     #select this
TSR = 8 # [-]     #select this too
CP = 0.46     #   select assumed Cp
cut_in = 3 # [m/s]
cut_off = 25 # [m/s]
#////////////////////////////////

# Calculate radius of multi rotor
#radius = equivalent_radius / (n_rotors ** 0.5)
AREA = P_RATED/(CP*0.5*rho*V_RATED**3)
radius = np.sqrt(AREA/(np.pi*n_rotors))

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

def find_all_x_crossings(x_data, y_data, y_targets):
    crossings = {y: [] for y in y_targets}

    # Iterate over each segment in the data
    for i in range(len(x_data) - 1):
        x0, x1 = x_data[i], x_data[i + 1]
        y0, y1 = y_data[i], y_data[i + 1]
        
        # Check each target y-value
        for y in y_targets:
            # Check if the y-value crosses between y0 and y1
            if (y0 <= y <= y1) or (y1 <= y <= y0):
                # Linear interpolation to find the x value of the crossing
                if y0 != y1:  # Avoid division by zero
                    x_cross = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
                    crossings[y].append(x_cross)
    
    return crossings

# Example usage
a = np.arange(-.5,1,.01)
x_data = a
# y_data = CTglauert*(1-a)
y_data = 4*a*(1-a)**2
y_targets = [CP*1.2]

crossings = find_all_x_crossings(x_data, y_data, y_targets)
estimated_a = np.min(crossings[y_targets[0]])

CT = 4*estimated_a*(1-estimated_a)
T_RATED = CT*0.5*rho*V_RATED**2*AREA
Q_RATED_perrotor = 0.5 * rho * CP * np.pi * radius ** 5 /(TSR ** 3) * (V_RATED * TSR / radius) ** 2
T_RATED_perrotor = T_RATED/n_rotors
P_RATED_perrotor = P_RATED/n_rotors
print(f'{P_RATED=}, {T_RATED=}, {Q_RATED_perrotor=}, {T_RATED_perrotor=}, {P_RATED_perrotor=}, {radius=}, {CT=}, {CP=}')
print(f'{estimated_a=}')
