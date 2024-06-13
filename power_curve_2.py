'''
Improved way of power curve
'''

import numpy as np
from matplotlib import pyplot as plt
import rewrittensimpleBEMmodel as BEM
import inputs_BEM_powerCurve as inps
# +----------------------------+
# |          Inputs            |
# +----------------------------+
#equivalent_radius = 170 # [m]
rho = inps.rho # [kg/m^3]
n_rotors = inps.n_rotors
P_RATED = inps.P_RATED # [W]
V_RATED = inps.V_RATED # [m/s]     #select this
TSR = inps.TSR # [-]     #select this too
CP = BEM.CP
if np.isnan(CP):
    CP = inps.assumed_CP #   select assumed Cp
cut_in = inps.cut_in # [m/s]
cut_off = inps.cut_off # [m/s]
#////////////////////////////////

# Calculate radius of multi rotor
#radius = equivalent_radius / (n_rotors ** 0.5)

radius = BEM.Radius
if np.isnan(BEM.CP):
    radius = inps.Radiuss
AREA = n_rotors*np.pi*radius**2
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
estimated_a = BEM.a_total

CT = BEM.CT
if np.isnan(CT):
    estimated_a = np.min(crossings[y_targets[0]])
    CT = 4*estimated_a*(1-estimated_a)

# Initialize Array
U_array = np.linspace(0, 27, 350)

# Calculate Power and Torque
P_array = 0.5 * rho * CP * U_array ** 3 * np.pi * radius ** 2 * n_rotors
Q_array = 0.5 * rho * CP * np.pi * radius ** 5 /(TSR ** 3) * (U_array * TSR / radius) ** 2
CP_array = np.ones(len(U_array))*CP
# T_array = 0.5*rho*CT*U_array**2*np.pi * radius**2*n_rotors
# Apply cut-in, rated and cut-off constraints
P_array[U_array < cut_in] = 0
Q_array[U_array < cut_in] = 0
CP_array[U_array<cut_in] = 0
# T_array[U_array<cut_in] = 0
Q_array[P_array > P_RATED] = Q_array[P_array > P_RATED][0]
CP_array[P_array>P_RATED] = P_RATED/(0.5*rho*U_array[P_array>P_RATED]**3*np.pi*radius**2*n_rotors)
CP1_array = CP_array
# T_array[P_array > P_RATED] = P_RATED / U_array
P_array[P_array > P_RATED] = P_RATED
P_array[U_array > cut_off] = 0
Q_array[U_array > cut_off] = 0
CP_array[U_array > cut_off] = 0


P_array /= 1e6
Q_array /= 1e3
fig, axs = plt.subplots(1, 2)
fig.set_figheight(3)
fig.set_figwidth(10)
if __name__ == "__main__":
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


    plt.plot(U_array, CP_array, label='CP with varying speeds', color='red')
    plt.legend()
    plt.xlabel('Uinf')
    plt.ylabel('CP')
    plt.title('CP with varying speeds')
    plt.show()
import pickle
with open('interpolated_function.pkl', 'rb') as f:
    loaded_interp_function = pickle.load(f)

CT_array = []
# from scipy.optimize import fsolve
# for ix in range(len(U_array)):
#     desired_CP = CP_array[ix]
#     def equation(x):
#         return loaded_interp_function(x) - desired_CP
#     # Use fsolve to find the root
#     initial_guess = inps.pitch  # Initial guess for the root
#     pitch_necessary = fsolve(equation, initial_guess)
#     CP_graph, CT_graph = BEM.BEMsolver_ale(pitch_necessary)
#     CT_array.append(CT_graph)

# for ix in range(len(CP1_array)):
#     x_data1 = np.arange(-10,25,1)
#     # y_data = CTglauert*(1-a)
#     y_data1 = loaded_interp_function(x_data1)
#     y_targets1 = CP1_array[ix]

#     crossings1 = find_all_x_crossings(x_data1, y_data1, y_targets1)
#     print(crossings1)
#     # estimated_pitch = np.min(crossings1[y_targets1[0]])
#     estimated_pitch = crossings1[y_targets1[0]]
#     print(estimated_pitch)
#     CP_graph, CT_graph = BEM.BEMsolver_ale(estimated_pitch)
#     CT_array.append(CT_graph)
# CT_array = np.array(CT_array)
# plt.plot(U_array, CT_array, label='CT with varying speeds', color='red')
# plt.legend()
# plt.xlabel('Uinf')
# plt.ylabel('CT')
# plt.title('CT with varying speeds')
# plt.show()

T_RATED = CT*0.5*rho*V_RATED**2*AREA
Q_RATED_perrotor = 0.5 * rho * CP * np.pi * radius ** 5 /(TSR ** 3) * (V_RATED * TSR / radius) ** 2
T_RATED_perrotor = T_RATED/n_rotors
P_RATED_perrotor = P_RATED/n_rotors
print(f'{P_RATED=}, {T_RATED=}, {Q_RATED_perrotor=}, {T_RATED_perrotor=}, {P_RATED_perrotor=}, {radius=}, {CT=}, {CP=}')
print(f'{estimated_a=}')
CP_glauert = 4*estimated_a*(1-estimated_a)**2
print((CP_glauert-CP)/CP_glauert)
print("speed at rotors is: ", V_RATED*(1-estimated_a))
diameter = radius*2
print(f'{diameter=}')
print(CP*0.5*rho*n_rotors*np.pi*radius**2*V_RATED**3)



# Uinfinity = np.arange(10, )
# BEM.BEMsolver_mark()