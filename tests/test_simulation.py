import pytest
import numpy as np


import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from horseshoe_vortex_model.simulation import Simulation, Horseshoe


@pytest.fixture(name='single_horseshoe')
def simple_horseshoe_sim():
    cor_1 = np.array((
        1, 0, 0
    ))
    cor_2 = np.array((
        1, 2, 0
    ))
    h_s = Horseshoe.initialize_horseshoe(cor_1, cor_2, 1)
    sim = Simulation([h_s], 1.0)
    return sim


def test_deflection_in_plane(single_horseshoe) -> None:
    pos = single_horseshoe.simulate([np.array((2, 1, 0))], 1, 1)
    res = pos[0]

    expected = np.array(
        (
            3,
            1,
            1/(4*np.pi) * (4/(2**0.5) + 2)
        )
    )

    assert all(res[i] == expected[i] for i in range(3))
