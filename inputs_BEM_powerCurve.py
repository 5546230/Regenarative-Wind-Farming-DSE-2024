import numpy as np

V_RATED = 13.48  #13.6 is effectively the biggest for 33 rotors
TSR = 8
P_RATED = 30*10**6
n_rotors = 33
rho = 1.225
cut_in = 3 # [m/s]
cut_off = 25 # [m/s]
pitch = 2 #degrees  #tip pitch i think
iteration = True
init_Radius = 31
tip_chord = 1
root_chord = 3
root_twist = -14
optimize = False
assumed_CP = 0.46
automatic_radius = True


AREA = P_RATED/(assumed_CP*0.5*rho*V_RATED**3)
#print(f'{AREA=}')
Radiuss = np.sqrt(AREA/(np.pi*n_rotors))
print(f'{Radiuss=}')

ale_shit = False

if automatic_radius:
    init_Radius = Radiuss

