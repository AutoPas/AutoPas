"""A hot periodic Lennard-Jones crystal melting, with the same physics computed natively and in Python.

Run with: python lennardJonesCube.py
"""

import time

import matplotlib.pyplot as plt
import numpy as np

import autopas as ap
from autopas.plotting import plot_history, plot_particles

TEMPERATURE = 2.0


class PythonLennardJones(ap.PairwiseFunctor):
    """12-6 Lennard-Jones, truncated at the cutoff, for sigma = epsilon = 1. Does what ap.LennardJones does."""

    def __init__(self, cutoff):
        super().__init__()
        self.shift = 4.0 * (cutoff**-12 - cutoff**-6)
        self.energy = 0.0

    def init_traversal(self):
        self.energy = 0.0

    def aos(self, i, j, newton3):
        dr = i.r - j.r
        inv_r2 = 1.0 / dr.dot(dr)
        inv_r6 = inv_r2**3
        force = 24.0 * inv_r2 * inv_r6 * (2.0 * inv_r6 - 1.0) * dr
        i.add_f(force)
        if newton3:
            j.sub_f(force)
        # Pairs reaching across the boundary are seen by their owned partner only. They count half.
        weight = (i.is_owned + (j.is_owned if newton3 else 0)) / 2.0
        self.energy += weight * (4.0 * inv_r6 * (inv_r6 - 1.0) - self.shift)

    def potential_energy(self):
        return self.energy


def make_simulation(functor, iterations):
    """The configuration is a dictionary with the same keys as md-flexible's YAML files."""
    sim = ap.AutoPas(
        {
            "container": ["LinkedCells"],
            "traversal": ["lc_c08"],
            "data-layout": ["AoS", "SoA"],
            "newton3": ["enabled"],
            "cutoff": 2.5,
            "verlet-skin-radius": 0.3,
            "verlet-rebuild-frequency": 10,
            "deltaT": 0.002,
            "iterations": iterations,
            "boundary-type": ["periodic"] * 3,
            "box-min": [0.0, 0.0, 0.0],
            "box-max": [9.6, 9.6, 9.6],
            "log-level": "warn",
            "Sites": {0: {"mass": 1.0, "epsilon": 1.0, "sigma": 1.0}},
        },
        functors=functor,
    )
    # A generator yields plain arrays, so any initial state can be derived from it. Here: thermal velocities.
    crystal = ap.generators.CubeGrid([8, 8, 8], particle_spacing=1.2, bottom_left_corner=[0.1, 0.1, 0.1]).generate()
    velocities = np.random.default_rng(42).normal(0.0, np.sqrt(TEMPERATURE), crystal["positions"].shape)
    sim.add_particles(crystal["positions"], velocities=velocities)
    return sim


def main():
    # LennardJones and not LennardJonesHWY, because the plot of the energies needs correct values for all data layouts.
    native = make_simulation(ap.LennardJones(), iterations=1000)
    start = time.perf_counter()
    native()  # Runs the iterations from the configuration.
    native_time = time.perf_counter() - start
    native(500)  # Resumes from the current state.

    python = make_simulation(PythonLennardJones(2.5), iterations=20)
    start = time.perf_counter()
    python()
    python_time = time.perf_counter() - start

    print(native)
    print(f"Native: {native_time / 1000 * 1e3:8.2f} ms per iteration, Python: {python_time / 20 * 1e3:8.2f} ms")
    print(f"Temperature: {native.temperature:.3f}, Configuration: {native.current_configuration}")
    print(f"Forces agree: {np.allclose(make_simulation(ap.LennardJones(), 20)().forces, python.forces, atol=1e-8)}")

    plot_particles(native, color_by="speed")
    plot_history(native)
    plt.show()


if __name__ == "__main__":
    main()
