from Structures_3 import Sizing, Steel
import numpy as np


if __name__ == '__main__':
    MAT = Steel()
    L = 10
    t = 0.01
    R_out = 5

    P_axial = 1000
    P_bending = 1000
    l_bending = 5

    TEST = Sizing(length=L, material=MAT)


    'MASS'
    expected = 7800*np.pi
    returned = TEST.calc_mass(t=t, R_out=R_out)
    print(f'mass={returned:.2f}, ratio={returned/expected}')

    'AREA'
    expected = np.pi/10
    returned = TEST.calc_area(t=t, R_out=R_out)
    print(f'mass={returned:.2f}, ratio={returned/expected}')

    'INERTIA'
    expected = 5*np.pi / 4
    returned = TEST.calc_I(t=t, R_out=R_out)
    print(f'mass={returned:.2f}, ratio={returned / expected}')

    'AXIAL THICKNESS'
    expected = P_axial/(2*np.pi*R_out*MAT.sig_y)
    returned = TEST.axial_thickness(axial_force=P_axial, R_out=R_out)
    print(f't_axial={returned:.2f}, ratio={returned / expected}')

    'BENDING THICKNESS'
    expected = P_bending*l_bending/(np.pi*R_out**2*MAT.sig_y)
    returned = TEST.bending_thickness(point_force=P_bending, position=l_bending, R_out=R_out)
    print(f't_axial={returned:.2f}, ratio={returned / expected}')

    'BUCKLING THICKNESS'
    # WRONG

