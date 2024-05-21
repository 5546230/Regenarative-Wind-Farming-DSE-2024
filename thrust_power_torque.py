# this program will calculate thrust power and torque

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Variables to change
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\
magic_radius_number_providedbyTiago = 170
rho = 1.225 #SI units

nrotors = 33
radius = magic_radius_number_providedbyTiago/((nrotors*(1+0))**0.5) #si units
Uinf = 10#si units
Power = 30 *10**6 #Watts
TSR = 8
max_Cp = 0.5 # np.max(CTglauert*(1-a)
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def CTfunction(a, glauert = False):
    """
    This function calculates the thrust coefficient as a function of induction factor 'a'
    'glauert' defines if the Glauert correction for heavily loaded rotors should be used; default value is false
    """
    CT = np.zeros(np.shape(a))
    CT = 4*a*(1-a)  
    if glauert:
        CT1=1.816
        a1=1-np.sqrt(CT1)/2
        CT[a>a1] = CT1-4*(np.sqrt(CT1)-1)*(1-a[a>a1])
    
    return CT
  
    
def ainduction(CT):
    """
    This function calculates the induction factor 'a' as a function of thrust coefficient CT 
    including Glauert's correction
    """
    a = np.zeros(np.shape(CT))
    CT1=1.816
    CT2=2*np.sqrt(CT1)-CT1
    a[CT>=CT2] = 1 + (CT[CT>=CT2]-CT1)/(4*(np.sqrt(CT1)-1))
    a[CT<CT2] = 0.5-0.5*np.sqrt(1-CT[CT<CT2])
    return a
a = np.arange(-.5,1,.01)
CTmom = CTfunction(a) # CT without correction
CTglauert = CTfunction(a, True) # CT with Glauert's correction
a2 = ainduction(CTglauert)

fig1 = plt.figure(figsize=(12, 6))
plt.plot(a, CTmom, 'k-', label='$C_T$')
plt.plot(a, CTglauert, 'b--', label='$C_T$ Glauert')
plt.plot(a, CTglauert*(1-a), 'g--', label='$C_P$ Glauert')
plt.xlabel('a')
plt.ylabel(r'$C_T$ and $C_P$')
plt.grid()
plt.legend()
plt.show()

# print(np.max(CTglauert*(1-a)))


CP = Power/nrotors/(0.5*rho*Uinf**3*np.pi*radius**2)  #necessary CP.

print(CP)
if CP> max_Cp:
    CP = min(CP, max_Cp)
    Power = CP*(0.5*rho*Uinf**3*np.pi*radius**2)*nrotors

Uinf_graph = np.arange(0, 30, 0.01)
Power_graph_lst = []
# Torque_graph_lst = []
for Uinfinity in Uinf_graph:
    Power_graph = 30*10**6
    CP_graph = Power_graph/nrotors/(0.5*rho*Uinfinity**3*np.pi*radius**2)  #necessary CP.
    
    if CP_graph> max_Cp:
        CP_graph = min(CP, max_Cp)
        Power_graph = CP*(0.5*rho*Uinfinity**3*np.pi*radius**2)*nrotors
        
    if Uinfinity > 25 or Uinfinity <3:   #cut in and out speed
        Power_graph = 0
        Torque_graph = 0
    # omega_graph =  Uinfinity*TSR/radius
    # Torque_graph = Power_graph / omega_graph
    Power_graph_lst.append(Power_graph)
    # Torque_graph_lst.append(Torque_graph)

Power_graph_lst = np.array(Power_graph_lst)
# Torque_graph_lst = np.array(Torque_graph_lst)

fig2 = plt.figure(figsize=(12, 6))
plt.plot(Uinf_graph, Power_graph_lst)
plt.xlabel('Wind Speed [m/s]')
plt.ylabel('Power Generated [W]')
plt.grid()
plt.show()

# fig3 = plt.figure(figsize=(12, 6))
# plt.plot(Uinf_graph, Torque_graph_lst)
# plt.xlabel('Wind Speed')
# plt.ylabel('Torque Generated')
# plt.grid()
# plt.show()
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
x_data = a
# y_data = CTglauert*(1-a)
y_data = 4*a*(1-a)**2
y_targets = [CP]

crossings = find_all_x_crossings(x_data, y_data, y_targets)
for y, x_vals in crossings.items():
    print(f"The graph crosses y = {y} at x = {x_vals}")
    print(np.min(crossings[CP]))

necessary_a = a[np.argmin(np.abs(CTglauert*(1-a)-CP))]
print("necessary a is", necessary_a)
necessary_a = np.min(crossings[CP])   #necessary a to get the necessary CP. 
print(" updated necessary a is ", necessary_a)
Area = np.pi*radius**2*nrotors

thrustperrotor = Power/Uinf/(1-necessary_a)/nrotors
print("Thrust per rotor is ", thrustperrotor, " at induced speed of ", Uinf*(1-necessary_a))
CT = thrustperrotor/(0.5*rho*Uinf**2*np.pi*radius**2)
print("CT is ", CT, " which should be the same as ", 4*necessary_a*(1-necessary_a))

omega =  Uinf*TSR/radius
Torqueperrotor = Power / omega/nrotors
print("Torque per rotor is ", Torqueperrotor)
print("Power per rotor is ", Power/nrotors)
print("Total Power is: ", Power)
print("total thrust", thrustperrotor*nrotors)

print("Cp is ", CP, "which should be the same as", 4*necessary_a*(1-necessary_a)**2)

print("omega is ", omega)