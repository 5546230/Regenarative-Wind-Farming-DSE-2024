from collections.abc import Iterable

import cProfile

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize


def normalise(vec: np.ndarray) -> np.ndarray:
    '''
    Normalise a vector

    Args:
        vec (np.ndarray): The vector to be normalised

    Returns:
        np.ndarray: The normalised vector
    '''
    return vec/np.linalg.norm(vec, axis=1)[:, None]


VORTEX_ELEMENT_DTYPE = np.dtype({
    'names': ('cor_1', 'cor_2', 'strength'),
    'formats': ('(3,)f8', '(3,)f8', 'f8')
})


class Horseshoe:
    '''
    A class to represent a horseshoe vortex element

    Attributes:
        elements (VORTEX_ELEMENT_DTYPE): The elements of the horseshoe
    '''
    @staticmethod
    def initialize_horseshoe(corner_1: np.ndarray, corner_2: np.ndarray, strength: np.ndarray) -> VORTEX_ELEMENT_DTYPE:
        '''
        Initialize a horseshoe vortex element

        Args:
            corner_1 (np.ndarray): The first corner of the horseshoe
            corner_2 (np.ndarray): The second corner of the horseshoe
            strength (np.ndarray): The strength of the horseshoe

        Returns:
            np.ndarray: The horseshoe vortex element
        '''
        inf_1 = np.column_stack(
            (10e10 * np.ones(corner_1.shape[0]), corner_1[:, 1:]))
        inf_2 = np.column_stack(
            (10e10 * np.ones(corner_2.shape[0]), corner_2[:, 1:]))

        element_1 = np.array(
            list(zip(corner_1, inf_1, strength)), dtype=VORTEX_ELEMENT_DTYPE)
        element_2 = np.array(
            list(zip(corner_1, corner_2, -strength)), dtype=VORTEX_ELEMENT_DTYPE)
        element_3 = np.array(
            list(zip(corner_2, inf_2, -strength)), dtype=VORTEX_ELEMENT_DTYPE)

        elements = np.column_stack((element_1, element_2, element_3))

        horseshoes = np.array(elements.flatten(), dtype=VORTEX_ELEMENT_DTYPE)
        return horseshoes

    def __str__(self) -> str:
        return str(self.elements[1])


class Simulation:
    '''	
    A class to represent a simulation of a wind system	

    Attributes:
        horseshoes (VORTEX_ELEMENT_DTYPE): The horseshoes of the wind system
        freestream (float): The velocity of the freestream
        height (float): The height of the wind system
    '''

    def __init__(self, horseshoes: VORTEX_ELEMENT_DTYPE, velocity: float, height: float = 250) -> None:
        self.horseshoes = horseshoes
        self.freestream = velocity
        self.height = height

    def get_streamlines(self, initial_pos: np.ndarray, dt: float, t_tot: float = 1) -> solve_ivp:
        '''
        Get the streamlines of the wind system for a given initial position

        Args:
            initial_pos (np.ndarray): The initial position of the particle
            dt (float): The time step of the simulation
            t_tot (float): The total time of the simulation

        Returns:
            solve_ivp: The solution of the simulation
        '''
        horseshoes = self.horseshoes
        ab = horseshoes['cor_2'] - horseshoes['cor_1']
        ab_n = np.linalg.norm(ab, axis=1)
        ab_norm = normalise(ab)
        n_horseshoes = len(horseshoes)

        def velocity(t: float, pos: np.ndarray) -> np.ndarray:
            '''
            Calculate the velocity of a particle at a given position

            Args:
                t (float): The time
                pos (np.ndarray): The position of the particle

            Returns:
                np.ndarray: The velocity of the particle
            '''
            pos = np.repeat(pos, n_horseshoes, axis=1).T
            ac = pos - horseshoes['cor_1']
            bc = pos - horseshoes['cor_2']
            ac_norm = normalise(ac)
            bc_norm = normalise(bc)
            cross_ac_bc = np.cross(ac, bc)
            h = np.linalg.norm(cross_ac_bc, axis=1) / ab_n
            phi_a = np.clip(np.einsum('ij,ij->i', ac_norm, ab_norm), -1, 1)
            phi_b = np.clip(np.einsum('ij,ij->i', -ab_norm, bc_norm), -1, 1)
            v = horseshoes['strength'] / \
                (4 * np.pi * h) * (phi_a + phi_b)
            direction = np.cross(ab, ac)
            direction = direction / np.linalg.norm(direction, axis=1)[:, None]
            output_vec = np.sum(direction * v[:, None], axis=0)
            output_vec[0] += self.freestream
            return output_vec

        def reached_top(t: float, pos: np.ndarray,
                        height: float) -> float: return pos[2] - height

        def reached_top_event(t: float, pos: np.ndarray) -> float: return reached_top(
            t, pos, self.height)

        reached_top_event.terminal = True
        reached_top_event.direction = 1

        sol = solve_ivp(velocity, [0, t_tot], initial_pos, t_eval=np.arange(
            0, t_tot, dt), events=[reached_top_event], vectorized=True)

        return sol

    def simulate(self, initial_positions: Iterable[np.ndarray], dt: float, t_tot: float = 1) -> list[solve_ivp]:
        '''
        Simulate the wind system for a given set of initial positions

        Args:
            initial_positions (Iterable[np.ndarray]): The initial positions of the particles
            dt (float): The time step of the simulation
            t_tot (float): The total time of the simulation

        Returns:
            list[solve_ivp]: The solutions of the simulations
        '''
        sols = list(self.get_streamlines(pos, dt, t_tot)
                    for pos in initial_positions)

        return sols

    def plot_horseshoes(self) -> None:
        fig, ax = plt.figure()


class WindSystem:
    '''
    A class to represent a wind system

    Attributes:
        horseshoes (VORTEX_ELEMENT_DTYPE): The horseshoes of the wind system
        height (float): The height of the wind system
        velocity (float): The velocity of the wind system
        dt (float): The time step of the simulation
    '''

    def __init__(self, horseshoes: VORTEX_ELEMENT_DTYPE, height: float, velocity: float, dt: float) -> None:
        self.velocity = velocity
        self.height = height
        self.dt = dt
        self.horseshoes = horseshoes

    @property
    def horseshoes(self) -> Iterable[Horseshoe]:
        return self._horseshoes

    @horseshoes.setter
    def horseshoes(self, val: VORTEX_ELEMENT_DTYPE) -> None:
        '''
        Set the horseshoes of the wind system

        Args:
            val (VORTEX_ELEMENT_DTYPE): The horseshoes to be set
        '''
        self._horseshoes = val

        cor_1_below = self.horseshoes[1::3]['cor_1'].copy()
        cor_1_below[:, 2] *= -1

        cor_2_below = self.horseshoes[1::3]['cor_2'].copy()
        cor_2_below[:, 2] *= -1

        strengths = self.horseshoes[1::3]['strength']

        self._below_horseshoes = Horseshoe.initialize_horseshoe(
            cor_1_below, cor_2_below, strengths)

        self.total_horseshoes = np.concatenate(
            (self._horseshoes, self._below_horseshoes))

        self.simulation = Simulation(
            self.total_horseshoes, self.velocity, self.height)

    def find_circulation(self, distance: float) -> float:
        '''
        Find the circulation required to reach a certain distance from the initial position to the top of the domain

        Args:
            distance (float): The distance to be reached

        Returns:
            float: The circulation required to reach the distance
        '''

        def objective(circ: float) -> float:
            '''
            Objective function to minimize the difference between the distance found and the desired distance

            Args:
                circ (float): The circulation to be used

            Returns:
                float: The difference between the distance found and the desired distance
            '''
            outcome = self.find_distance_circ(circ)
            return (outcome - distance) ** 2

        initial = self.horseshoes[0]['strength']

        optimized = minimize(objective, initial)

        return optimized.x[0]

    def find_distance_circ(self, circulation: float) -> float:
        '''
        Find the distance from the initial position to the top of the domain with a given circulation

        Args:
            circulation (float): The circulation to be used

        Returns:
            float: The distance from the initial position to the top of the domain
        '''

        # Update the strengths
        self.total_horseshoes['strength'] *= circulation / \
            self.total_horseshoes['strength'][0]

        return self.find_distance()

    def find_distance(self) -> float:
        '''
        Find the distance from the initial position to the top of the domain
        '''

        initial_pos = np.array((
            0,
            100,
            10
        ))

        sol = self.simulation.get_streamlines(initial_pos, self.dt, 1500)

        return sol.y_events[0][0][0]


def main_1() -> None:

    # TODO: Check numbers (height, span)
    heights = np.array([235, 165, 95, 25, -235, -165, -95, -25])
    span = 250
    gamma = 270
    height = 250

    corner_1 = np.column_stack(
        (np.ones(heights.shape[0]), 10 * np.ones(heights.shape[0]), heights))
    corner_2 = np.column_stack(
        (np.ones(heights.shape[0]), (10 + span) * np.ones(heights.shape[0]), heights))
    strength = gamma * (-1 + 2 * (heights > 0))

    vortices = Horseshoe.initialize_horseshoe(corner_1, corner_2, strength)

    # Create a simulation instance
    simulation = Simulation(
        horseshoes=vortices, velocity=9.29)

    fig = plt.figure()
    three_d = True
    if three_d:
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.add_subplot(111)

    y_values = np.array([50, 100, 150])
    z_values = np.array([10, 50, 100, 150, 200])

    y, z = np.meshgrid(y_values, z_values)
    initial_positions = np.dstack((np.zeros_like(y), y, z)).reshape(-1, 3)

    sols = simulation.simulate(initial_positions, dt=0.1, t_tot=500)
    for i, sol in enumerate(sols):
        # print(sol.t_events, sol.y_events)
        pos = sol.y
        # if i == 0:
        # print(sol.y)

        if three_d:
            ax.plot(pos[0], pos[1],
                    pos[2], label='Streamline')
            ax.scatter(pos[0][0], pos[1][0], pos[2][0],
                       color='red', label='Starting point')
        else:
            ax.plot(pos[0][::10], pos[2][::10], label=f'Streamline {i}')

    # for i, vortex in enumerate(simulation.horseshoes):
    #     locs = np.array(vortex.locations)
    #     ax.plot(locs[:, 0], locs[:, 1], locs[:, 2], label=f'Vortex {i}')

    # print(x, y, z)

    # ax.quiver(
    #     x, y, z,
    #     simulation.velocity_mesh[:,:,:, 0],
    #     simulation.velocity_mesh[:,:,:, 1],
    #     simulation.velocity_mesh[:,:,:, 2],
    #     length = 0.5,
    #     normalize=True
    #     )

    # Set labels
    if three_d:
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
    else:
        ax.set_xlabel('X')
        ax.set_ylabel('Z')

    # ax.legend()
    plt.show()


def main_2() -> None:
    # Define a single vortex with one location
    # TODO: Check numbers (height, span)
    heights = np.array([235, 165, 95, 25])
    span = 250
    gamma = 270
    height = 250
    V = 9.29

    corner_1 = np.column_stack(
        (np.ones(heights.shape[0]), 10 * np.ones(heights.shape[0]), heights))
    corner_2 = np.column_stack(
        (np.ones(heights.shape[0]), (10 + span) * np.ones(heights.shape[0]), heights))
    strength = gamma * np.ones(heights.shape[0])

    vortices = Horseshoe.initialize_horseshoe(corner_1, corner_2, strength)

    system = WindSystem(vortices, height, V, 1)

    circ_req = system.find_circulation(2000)
    lift = circ_req * 1.225 * V * 250
    cl_c = lift/(0.5*1.225*V**2*250)
    print(f'{circ_req=}, {lift=}, {cl_c=}')
    print(system.find_distance_circ(250))


def analyze_profile(output_file: str) -> None:
    import pstats
    p = pstats.Stats(output_file)
    p.sort_stats('time').print_stats(10)


if __name__ == "__main__":
    # cProfile.run('main_2()', 'output.prof')
    # analyze_profile('output.prof')
    main_2()
