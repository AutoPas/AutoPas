"""Shared configurations for the tests of the Python interface."""

import matplotlib
import pytest

import autopas as ap

matplotlib.use("Agg")

CUTOFF = 2.5


def small_config(iterations=20, **overrides):
    """A periodic crystal with a fixed algorithm configuration, so that no tuning takes place."""
    config = {
        "container": ["LinkedCells"],
        "traversal": ["lc_c08"],
        "data-layout": ["AoS"],
        "newton3": ["enabled"],
        "cutoff": CUTOFF,
        "verlet-skin-radius": 0.3,
        "verlet-rebuild-frequency": 10,
        "deltaT": 0.002,
        "iterations": iterations,
        "boundary-type": ["periodic"] * 3,
        "box-min": [0.0, 0.0, 0.0],
        "box-max": [7.2, 7.2, 7.2],
        "log-level": "warn",
        "Sites": {0: {"mass": 1.0, "epsilon": 1.0, "sigma": 1.0}},
        "Objects": {
            "CubeGrid": {
                0: {
                    "particles-per-dimension": [6, 6, 6],
                    "particle-spacing": 1.2,
                    "bottomLeftCorner": [0.1, 0.1, 0.1],
                    "velocity": [0.0, 0.0, 0.0],
                    "particle-type-id": 0,
                }
            }
        },
    }
    config.update(overrides)
    return config


class PythonLennardJones(ap.PairwiseFunctor):
    """The physics of ap.LennardJones for sigma = epsilon = 1 and a truncated, shifted potential."""

    def __init__(self, cutoff=CUTOFF, **kwargs):
        super().__init__(**kwargs)
        self.shift = 4.0 * (cutoff**-12 - cutoff**-6)
        self.energy = 0.0
        self.calls = 0

    def init_traversal(self):
        self.energy = 0.0

    def aos(self, i, j, newton3):
        self.calls += 1
        dr = i.r - j.r
        inv_r2 = 1.0 / dr.dot(dr)
        inv_r6 = inv_r2**3
        force = 24.0 * inv_r2 * inv_r6 * (2.0 * inv_r6 - 1.0) * dr
        i.add_f(force)
        if newton3:
            j.sub_f(force)
        weight = (i.is_owned + (j.is_owned if newton3 else 0)) / 2.0
        self.energy += weight * (4.0 * inv_r6 * (inv_r6 - 1.0) - self.shift)

    def potential_energy(self):
        return self.energy


@pytest.fixture
def config():
    return small_config()
