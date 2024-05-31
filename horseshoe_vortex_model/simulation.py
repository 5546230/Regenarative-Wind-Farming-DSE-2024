from collections.abc import Iterable

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize


def normalise(vec):
    return vec/np.linalg.norm(vec)


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
    def __init__(self, elements: Iterable[VortexElement]) -> None:
        self.elements = elements

    @staticmethod
    def initialize_horseshoe(corner_1, corner_2, strength):
        inf_1 = np.array((10e10, *corner_1[1:]))
        inf_2 = np.array((10e10, *corner_2[1:]))

        element_1 = VortexElement(
            corner_1, inf_1, strength)
        element_2 = VortexElement(
            corner_1, corner_2, -strength)
        element_3 = VortexElement(
            corner_2, inf_2, -strength)

        horseshoe = Horseshoe([element_1, element_2, element_3])
        return horseshoe


class Simulation:
    def __init__(self, horseshoes: Iterable[Horseshoe], velocity, height=400) -> None:
        self.horseshoes = horseshoes
        self.freestream = velocity
        self.height = height

    def get_streamlines(self, initial_pos, dt, t_tot=1):
        def velocity(t, pos):
            output_vec = np.zeros((3, ))
            output_vec[0] = self.freestream
            for i, vortex in enumerate(self.horseshoes):
                for j, element in enumerate(vortex.elements):
                    # print(element, j)
                    ab = element[1] - element[0]
                    ac = pos - element[0]
                    bc = pos - element[1]

                    ab_norm = normalise(ab)
                    ac_norm = normalise(ac)
                    bc_norm = normalise(bc)

                    cross = np.linalg.norm(np.cross(ab, ac))
                    if np.isnan(cross):
                        print('nan')
                        print(element, pos, cross)
                        return
                    if np.isclose(cross, 0.):
                        # print(element)
                        continue

                    h = np.linalg.norm(np.cross(ac, bc))/np.linalg.norm(ab)
                    if j == 0:
                        v = element.strength / \
                            (4 * np.pi * h) * (np.dot(ac, ab) /
                                               (np.linalg.norm(ab) * np.linalg.norm(ac))+1)
                        # print(v)
                    elif j == 1:
                        phi_a = np.arccos(
                            np.dot(ac, ab) / (np.linalg.norm(ab) * np.linalg.norm(ac)))
                        phi_b = np.arccos(np.clip(
                            np.dot(-ab_norm, bc_norm), -1, 1))
                        if np.isnan(phi_b):
                            # print(ac, bc, np.dot(-ab_norm, bc_norm),
                            #       np.linalg.norm(ab), np.linalg.norm(bc))
                            return (0, 0, 0)
                        v = element.strength / \
                            (4 * np.pi * h) * (np.cos(phi_a) + np.cos(phi_b))
                        # print(v)

                    elif j == 2:
                        v = element.strength / \
                            (4 * np.pi * h) * \
                            (np.clip(np.dot(ab_norm, ac_norm), -1, 1)+1)
                        # print(v)

                    direction = np.cross(ab, ac)
                    # if np.isnan(np.linalg.norm(direction)):
                    # print('nan1')
                    # print(element, pos, direction)
                    direction = direction / np.linalg.norm(direction)

                    # if
                    output_vec += direction.reshape((3,)) * v
                    # print(output_vec, j)

                    # print(i, j)
            # print(output_vec)
            return output_vec

        def reached_top(t, pos, height):
            return pos[2] - height

        def reached_top_event(t, pos): return reached_top(t, pos, self.height)

        reached_top_event.terminal = True
        reached_top_event.direction = 1

        sol = solve_ivp(
            velocity, [0, t_tot], initial_pos, t_eval=np.arange(0, t_tot, dt), events=[reached_top_event])

        # return initial_pos + velocity(1, initial_pos)

        return sol

    def simulate(self, initial_positions, dt, t_tot=1):
        sols = list(self.get_streamlines(pos, dt, t_tot)
                    for pos in initial_positions)

        return sols


class WindSystem:
    def __init__(self, horseshoes: Iterable[Horseshoe], height: float, velocity: float, dt: float) -> None:
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
        self.simulation = Simulation(
            self.horseshoes, self.velocity, self.height)

    def find_circulation(self, distance: float):

        def objective(circ: float) -> float:
            outcome = self.find_distance_circ(circ)
            return (outcome - distance) ** 2

        initial = self.horseshoes[0].elements[0].strength

        optimized = minimize(objective, initial)

        return optimized.x[0]

    def find_distance_circ(self, circulation: float) -> float:
        self.horseshoes = list(
            Horseshoe.initialize_horseshoe(hs.elements[1].cor_1, hs.elements[1].cor_2, circulation) for hs in self.horseshoes
        )
        return self.find_distance()

    def find_distance(self) -> float:

        initial_pos = np.array((
            0,
            150,
            0
        ))

        sol = self.simulation.get_streamlines(initial_pos, self.dt, 1500)
        return sol.y_events[0][0][0]


def main_1():
    # Define a single vortex with one location
    vortex1 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 310.)), np.array((1., 310., 310.)), 278.7)
    vortex2 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 265.)), np.array((1., 310., 265.)), 278.7)
    vortex3 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 170.)), np.array((1., 310., 170.)), 278.7)
    vortex4 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 125.)), np.array((1., 310., 125.)), 278.7)
    vortices = [vortex1, vortex2, vortex3, vortex4]

    # Create a simulation instance
    simulation = Simulation(
        horseshoes=vortices, velocity=9.29)

    fig = plt.figure()
    three_d = False
    if three_d:
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.add_subplot(111)

    initial_positions = [
        np.array((0, 50, 0)),
        np.array((0, 150, 0)),
        np.array((0, 250, 0)),
        # np.array((0, 50, 50)),
        np.array((0, 150, 50)),
        # np.array((0, 250, 50)),
        # np.array((0, 50, 100)),
        np.array((0, 150, 100)),
        # np.array((0, 250, 10)),
        # np.array((0, 50, 150)),
        np.array((0, 150, 150)),
        # np.array((0, 250, 150)),
        # np.array((0, 50, 200)),
        np.array((0, 150, 200)),
        # np.array((0, 250, 200)),
        # np.array((0, 50, 250)),
        np.array((0, 150, 250)),
        # np.array((0, 250, 250)),
    ]

    sols = simulation.simulate(initial_positions, dt=0.0001, t_tot=1500)
    for i, sol in enumerate(sols):
        print(sol.t_events, sol.y_events)
        pos = sol.y

        if three_d:
            ax.plot(pos[0][::10], pos[1][::10],
                    pos[2][::10], label='Streamline')
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
    vortex1 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 310.)), np.array((1., 310., 310.)), 366)
    vortex2 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 265.)), np.array((1., 310., 265.)), 366)
    vortex3 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 170.)), np.array((1., 310., 170.)), 366)
    vortex4 = Horseshoe.initialize_horseshoe(
        np.array((1., 10., 125.)), np.array((1., 310., 125.)), 366)
    vortices = [vortex1, vortex2, vortex3, vortex4]

    system = WindSystem(vortices, 400, 9.29, 0.001)

    print(system.find_distance())
    print(system.find_distance_circ(278.7))

    print(system.find_circulation(1700))


if __name__ == "__main__":
    main_2()
