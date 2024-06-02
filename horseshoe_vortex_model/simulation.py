from collections.abc import Iterable

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize


def normalise(vec):
    return vec/np.linalg.norm(vec, axis=1)[:, None]


VORTEX_ELEMENT_DTYPE = np.dtype({
    'names': ('cor_1', 'cor_2', 'strength'),
    'formats': ('(3,)f8', '(3,)f8', 'f8')
})

HORSESHOE_DTYPE = np.dtype([
    ('elements', VORTEX_ELEMENT_DTYPE, (1, 3))
])


class VortexElement:
    def __init__(self, cor_1, cor_2, strength) -> None:
        self.cor_1 = cor_1
        self.cor_2 = cor_2
        self.strength = strength

    def __getitem__(self, key):
        return [self.cor_1, self.cor_2][key]

    def __str__(self) -> str:
        return f'VortexElement; {self.cor_1}, {self.cor_2}, {self.strength=}'


class Horseshoe:
    @staticmethod
    def initialize_horseshoe(corner_1, corner_2, strength):
        inf_1 = np.column_stack(
            (10e10 * np.ones(corner_1.shape[0]), corner_1[:, 1:]))
        inf_2 = np.column_stack(
            (10e10 * np.ones(corner_2.shape[0]), corner_2[:, 1:]))
        # print(corner_1, corner_2, inf_1, inf_2, strength)

        element_1 = np.array(
            list(zip(corner_1, inf_1, strength)), dtype=VORTEX_ELEMENT_DTYPE)
        element_2 = np.array(
            list(zip(corner_1, corner_2, -strength)), dtype=VORTEX_ELEMENT_DTYPE)
        element_3 = np.array(
            list(zip(corner_2, inf_2, -strength)), dtype=VORTEX_ELEMENT_DTYPE)
        # print(element_1, element_2, element_3)

        elements = np.column_stack((element_1, element_2, element_3))

        horseshoes = np.array(elements.flatten(), dtype=VORTEX_ELEMENT_DTYPE)
        return horseshoes

    def __str__(self) -> str:
        return str(self.elements[1])


class Simulation:
    def __init__(self, horseshoes, velocity, height=250) -> None:
        self.horseshoes = horseshoes
        self.freestream = velocity
        self.height = height

    def get_streamlines(self, initial_pos, dt, t_tot=1):
        def velocity(t, pos):
            # Convert list to numpy array for vectorized operations
            horseshoes = self.horseshoes

            pos = np.repeat(pos, len(self.horseshoes), axis=1).T

            # Calculate ab, ac, bc
            ab = horseshoes['cor_2'] - horseshoes['cor_1']
            ac = pos - horseshoes['cor_1']
            bc = pos - horseshoes['cor_2']

            # Normalize ab, ac, bc
            ab_norm = normalise(ab)
            ac_norm = normalise(ac)
            bc_norm = normalise(bc)

            # Calculate h
            h = np.linalg.norm(np.cross(ac, bc), axis=1) / \
                np.linalg.norm(ab, axis=1)

            # Calculate phi_a, phi_b for j%3 == 1
            phi_a_mask = np.arange(len(horseshoes)) % 3 == 1
            phi_a = np.arccos(
                np.clip(np.einsum('ij,ij->i', ac_norm[phi_a_mask], ab_norm[phi_a_mask]), -1, 1))
            phi_b = np.arccos(np.clip(
                np.einsum('ij,ij->i', -ab_norm[phi_a_mask], bc_norm[phi_a_mask]), -1, 1))

            # Calculate v for j%3 == 1
            v = np.zeros_like(h)
            v[phi_a_mask] = horseshoes[phi_a_mask]['strength'] / \
                (4 * np.pi * h[phi_a_mask]) * (np.cos(phi_a) + np.cos(phi_b))

            # Calculate v for j%3 == 2 or j%3 == 0
            v_mask = np.logical_or(np.arange(len(horseshoes)) %
                                   3 == 2, np.arange(len(horseshoes)) % 3 == 0)
            v[v_mask] = horseshoes[v_mask]['strength'] / (4 * np.pi * h[v_mask]) * (
                np.clip(np.einsum('ij,ij->i', ab_norm[v_mask], ac_norm[v_mask]), -1, 1) + 1)

            # Calculate direction
            direction = np.cross(ab, ac)
            direction = direction / np.linalg.norm(direction, axis=1)[:, None]

            # Calculate output_vec
            output_vec = np.sum(direction * v[:, None], axis=0)
            output_vec[0] += self.freestream

            return output_vec

        def reached_top(t, pos, height): return pos[2] - height

        def reached_top_event(t, pos): return reached_top(t, pos, self.height)

        reached_top_event.terminal = True
        reached_top_event.direction = 1

        sol = solve_ivp(velocity, [0, t_tot], initial_pos, t_eval=np.arange(
            0, t_tot, dt), events=[reached_top_event], vectorized=True)

        return sol

    def simulate(self, initial_positions, dt, t_tot=1):
        sols = list(self.get_streamlines(pos, dt, t_tot)
                    for pos in initial_positions)

        return sols

    def plot_horseshoes(self) -> None:
        fig, ax = plt.figure()


class WindSystem:
    def __init__(self, horseshoes, height: float, velocity: float, dt: float) -> None:
        self.velocity = velocity
        self.height = height
        self.dt = dt
        self.horseshoes = horseshoes

    @property
    def horseshoes(self) -> Iterable[Horseshoe]:
        return self._horseshoes

    @horseshoes.setter
    def horseshoes(self, val) -> None:
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

    def find_circulation(self, distance: float):

        def objective(circ: float) -> float:
            outcome = self.find_distance_circ(circ)
            return (outcome - distance) ** 2

        initial = self.horseshoes[0]['strength']

        optimized = minimize(objective, initial)

        # Check if the optimization was successful
        # if not optimized.success:
        #     raise Exception(
        #         f'Optimization failed: {optimized.message}')

        return optimized.x[0]

    def find_distance_circ(self, circulation: float) -> float:

        # Update the strengths
        self.total_horseshoes['strength'] *= circulation / \
            self.total_horseshoes['strength'][0]

        return self.find_distance()

    def find_distance(self) -> float:

        initial_pos = np.array((
            0,
            100,
            10
        ))

        sol = self.simulation.get_streamlines(initial_pos, self.dt, 1500)

        # Check if the streamline reached the top
        # if not sol.status:
        #     raise Exception(
        #         f'Streamline did not reach the top: {sol.message}')

        return sol.y_events[0][0][0]


def main_1():

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

    y_values = np.array([-10, 50, 100, 150])
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


def main_2():
    # Define a single vortex with one location
    # TODO: Check numbers (height, span)
    heights = np.array([235, 165, 95, 25])
    span = 250
    gamma = 270
    height = 250

    corner_1 = np.column_stack(
        (np.ones(heights.shape[0]), 10 * np.ones(heights.shape[0]), heights))
    corner_2 = np.column_stack(
        (np.ones(heights.shape[0]), (10 + span) * np.ones(heights.shape[0]), heights))
    strength = gamma * np.ones(heights.shape[0])

    vortices = Horseshoe.initialize_horseshoe(corner_1, corner_2, strength)

    V = 9.29
    system = WindSystem(vortices, height, V, 1)

    # print(system.find_distance())
    # print(system.find_distance_circ(278.7))
    try:
        circ_req = system.find_circulation(1700)
        lift = circ_req * 1.225 * V * 250
        cl_c = lift/(0.5*1.225*V**2*250)
        print(f'{circ_req=}, {lift=}, {cl_c=}')
    except Exception as e:
        print(e)
    print(system.find_distance_circ(250))


if __name__ == "__main__":
    main_2()
