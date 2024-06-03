from Truss.Honeycomb_Truss import Hexagonal_Truss
from drag_calc import Drag
import numpy as np


def case_1(hex: Hexagonal_Truss ,):

    'Thrust forces:'
    T_rated_per_rot = 119.5e3  # [N]
    front_T_indices = hex.find_midpoint_indices(side='front')
    back_T_indices = hex.find_midpoint_indices(side='back')

    T_force = np.array([[0] ,[T_rated_per_rot] ,[0]] ) /2
    front_T_forces = np.repeat(T_force, front_T_indices.size, axis=1)
    back_T_forces = np.repeat(T_force, back_T_indices.size, axis=1)
    # print(front_T_forces.shape, back_T_forces.shape)

    'Lift and Drag forces:'
    L_per_wing = 3e6  # [N]
    L_application = np.array([[0] ,[0] ,[-L_per_wing]] ) /4


    'Drag forces ...'
    drag = Drag()


    'Inertial forces ...'



    load_indices = np.hstack((front_T_indices, back_T_indices))
    loads = np.hstack((front_T_forces, back_T_forces))
    hex.load_indices = load_indices
    hex.applied_loads = loads
    return hex


if __name__ == "__main__":
    print('running analysis')