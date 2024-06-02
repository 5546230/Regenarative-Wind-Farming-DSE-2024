import pytest
import numpy as np


import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from horseshoe_vortex_model.simulation import Simulation, Horseshoe, WindSystem


@pytest.fixture(name='single_horseshoe')
def simple_horseshoe_sim():
    cor_1 = np.array((
        1, 0, 1
    ))
    cor_2 = np.array((
        1, 2, 1
    ))
    h_s = Horseshoe.initialize_horseshoe(cor_1, cor_2, 1)
    return h_s

@pytest.mark.skip
def test_deflection_in_plane(single_horseshoe) -> None:
    sim = Simulation([single_horseshoe], 1.0)
    pos = sim.simulate([np.array((2, 1, 1))], 1, 1)
    res = pos[0]
    print(res)

    expected = np.array(
        (
            3,
            1,
            1/(4*np.pi) * (4/(2**0.5) + 2) + 1
        )
    )
    print(expected)

    assert all(res[i] == expected[i] for i in range(3))


def test_deflection_xyplane():
    cor_1 = np.array((
        1, 0, 1
    ))
    cor_2 = np.array((
        1, 2, 1
    ))
    h_s = Horseshoe.initialize_horseshoe(cor_1, cor_2, 1)
    sys = WindSystem([h_s], 250, 1, 1)

    pos = sys.simulation.simulate([np.array((2, 1, 0))], 1, 1)
    res = pos[0]
    print(res)
