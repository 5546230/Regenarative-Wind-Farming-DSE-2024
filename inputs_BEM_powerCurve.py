import numpy as np

V_RATED = 10.59  #13.6 is effectively the biggest for 33 rotors
TSR = 8
P_RATED = 30*10**6
n_rotors = 33
rho = 1.225
cut_in = 3 # [m/s]
cut_off = 25 # [m/s]
pitch = -3.275370556176349#2 #degrees  #tip pitch i think
iteration = True     #iterates to converge for the correct radius. 
init_Radius = 19.89 #initial radius estimate. If you don't trust it, set automatic_radius to True, which estimates initial radius based on assumed_CP
automatic_radius = True   
tip_chord = 1
root_chord = 2#3
root_twist = -8.420070544909368#-14
optimize = True    #optimize blade geometry (pitch, twist, chords)
assumed_CP = 0.46  
ale_shit = False   #something bc alessandro wanted CT(TSR, pitch)
ale_shit_2 = False
printReynolds = False

AREA = P_RATED/(assumed_CP*0.5*rho*V_RATED**3)
#print(f'{AREA=}')
Radiuss = np.sqrt(AREA/(np.pi*n_rotors))
print(f'{Radiuss=}')



if automatic_radius:
    init_Radius = Radiuss

