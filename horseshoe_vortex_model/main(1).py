import time
from dataclasses import dataclass
import toml

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np


@dataclass
class Vec2D:
    y: float
    z: float


@dataclass
class Vortex:
    position: Vec2D
    strength: float


def compute_single_velocity(g2p: float, pt: Vec2D, mesh: np.ndarray) -> np.ndarray:
    v_field = np.zeros_like(mesh, dtype=np.float64)
    r_y = mesh[:, :, 0] - pt.y
    r_z = mesh[:, :, 1] - pt.z
    dist = np.hypot(r_y, r_z)
    v_field[:, :, 0] = r_z / (dist ** 2)
    v_field[:, :, 1] = -r_y / (dist ** 2)
    v_field[np.isclose(dist, 0.0), :] = 0
    v_field *= g2p
    return v_field


def compute_flow_field(mesh: np.ndarray, vor_list: list[Vortex]) -> np.ndarray:
    v_field = np.zeros_like(mesh, dtype=np.float64)
    for i, vor in enumerate(vor_list):
        v_vor = compute_single_velocity(vor.strength / (2 * np.pi), vor.position, mesh)
        v_field += v_vor
    return v_field


def mesh_create(ny: int, nz: int, origin: Vec2D, dy: float, dz: float) -> np.ndarray:
    y_scale = np.arange(ny) * dy + origin.y
    z_scale = np.arange(nz) * dz + origin.z
    yy, zz = np.meshgrid(y_scale, z_scale)
    return np.dstack((yy, zz))


if __name__ == '__main__':
    cfg_values = None
    # with open("horseshoe_vortex_model/config.toml", "rb") as f_in:
    cfg_values = toml.load("horseshoe_vortex_model/config.toml")
    sim_params = cfg_values["SIMULATION_PARAMETERS"]
    DT = float(sim_params["DT"])
    MAX_TIME_STEPS = int(sim_params["MAX_TIME_STEPS"])
    DOMAIN_W = int(sim_params["DOMAIN_W"])
    DOMAIN_H = int(sim_params["DOMAIN_H"])
    N_WIDTH = int(sim_params["N_WIDTH"])
    N_HEIGHT = int(sim_params["N_HEIGHT"])
    ORIGIN = Vec2D(*sim_params["ORIGIN"])
    sim_contents = cfg_values["SIMULATION_CONTENTS"]
    TRACKED_POINT = Vec2D(*sim_contents["TRACKED_POINT"])
    vor_parameters = sim_contents["VORTEX_LIST"]
    VORTEX_LIST = []
    for i, v_vor in enumerate(vor_parameters):
        VORTEX_LIST.append(
            Vortex(Vec2D(v_vor['position'][0], v_vor['position'][1]), v_vor['strength'])
        )

    dw = DOMAIN_W / (N_WIDTH - 1)
    dh = DOMAIN_H / (N_HEIGHT - 1)
    mesh_space = mesh_create(N_WIDTH, N_HEIGHT, ORIGIN, dw, dh)

    mesh_velocity = compute_flow_field(mesh_space, VORTEX_LIST)
    #   Simulate position of tracked point until it exits the domain
    MIN_Y = np.min(mesh_space[:, :, 0]);
    MAX_Y = np.max(mesh_space[:, :, 0]);
    MIN_Z = np.min(mesh_space[:, :, 1]);
    MAX_Z = np.max(mesh_space[:, :, 1]);
    MIN_LIM = np.array([MIN_Y, MIN_Z]);
    MAX_LIM = np.array([MAX_Y, MAX_Z])

    def is_pt_in_domain(r: np.ndarray) -> bool:
        return np.all((MIN_LIM < r) * (r < MAX_LIM))


    pos = np.array([TRACKED_POINT.y, TRACKED_POINT.z], dtype=np.float64)
    n = 0
    p = np.zeros((MAX_TIME_STEPS, 2), dtype=np.float64)
    p[0, :] = pos
    N_VORTEX = len(VORTEX_LIST)
    pos_vortex = np.zeros((N_VORTEX, 2), dtype=np.float64)
    str_vortex = np.zeros((N_VORTEX, 2), dtype=np.float64)
    for i, vor in enumerate(VORTEX_LIST):
        pos_vortex[i, :] = (vor.position.y, vor.position.z)
        str_vortex[i, :] = vor.strength / (2 * np.pi)

    t0 = time.perf_counter_ns()

    while is_pt_in_domain(pos) and n < MAX_TIME_STEPS:
        pos_relative = pos_vortex - pos
        dist_relative = np.hypot(pos_relative[:, 0], pos_relative[:, 1]) ** 2
        valid = np.logical_not(np.isclose(dist_relative, 0))
        contribution = pos_relative[valid, :] * str_vortex[valid, :]
        v_vor = np.zeros((np.sum(valid), 2), dtype=np.float64)
        v_vor[:, 0] = -contribution[:, 1] / dist_relative[valid]
        v_vor[:, 1] = +contribution[:, 0] / dist_relative[valid]
        pos += DT * np.sum(v_vor, axis=0)
        n += 1
        p[n, :] = pos
    t1 = time.perf_counter_ns()

    # print(n, "calculations took", (t1 - t0) / 1e6, "ms")

    #   Displays the velocity flow field along with the position of the tracked point
    fig: plt.Figure
    ax: plt.Axes
    fig, ax = plt.subplots(1, 1)
    ax.contourf(mesh_space[:, :, 0], mesh_space[:, :, 1], np.hypot(mesh_velocity[:, :, 0], mesh_velocity[:, :, 1]), norm=matplotlib.colors.LogNorm())
    ax.streamplot(mesh_space[:, :, 0], mesh_space[:, :, 1], mesh_velocity[:, :, 0], mesh_velocity[:, :, 1])
    ax.plot(p[:n, 0], p[:n, 1], color="red")
    ax.set_title("Velocity profile cross-section")
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_aspect("equal")
    # ax.set_xlim(-270, 270)
    # ax.set_ylim(-310, 310)
    plt.show()
