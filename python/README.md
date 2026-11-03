# AutoPas for Python

Python interface of [AutoPas](https://github.com/AutoPas/AutoPas), built with pybind11 on top of the molecular
dynamics application library and md-flexible.
It covers single-site molecular dynamics on one MPI rank, parallelized with OpenMP.

## Installation

```shell
conda create -n autopas -c conda-forge python=3.12 numpy matplotlib pytest
conda activate autopas
pip install ".[plot]"
```

This compiles AutoPas including all of its traversals for each functor, which took about a minute on a 10 core laptop.
The module is compiled with `-march=native`, so it only runs on the machine (type) it was built on.
Options of the C++ build are passed through pip:

```shell
pip install . -C cmake.define.AUTOPAS_PYTHON_FUNCTOR_AVX=ON
```

On macOS, the OpenMP runtime has to be available, e.g. via `conda install -c conda-forge llvm-openmp`.

## Usage

The notebook [`fallingDrop.ipynb`](../notebooks/fallingDrop.ipynb) walks through a whole simulation, including a GIF of it.
[`solarSystem.ipynb`](../notebooks/solarSystem.ipynb) simulates the solar system with a gravity functor written in
Python.
They additionally need `jupyterlab`.

```python
import autopas as ap
from autopas.plotting import plot_history, plot_particles

sim = ap.AutoPas(
    {
        "cutoff": 2.5,
        "deltaT": 0.002,
        "iterations": 100,
        "box-min": [0, 0, 0],
        "box-max": [9.6, 9.6, 9.6],
        "Sites": {0: {"mass": 1.0, "epsilon": 1.0, "sigma": 1.0}},
        "Objects": [ap.generators.CubeGrid([8, 8, 8], particle_spacing=1.2, bottom_left_corner=[0.1, 0.1, 0.1])],
    },
    functors=ap.LennardJonesHWY(),
)

sim()      # Starts the simulation and runs the 100 iterations from the configuration.
sim(50)    # Resumes from the current state for 50 more iterations.
sim(50, callback=lambda s: print(s.iteration, s.temperature), callback_interval=10)

sim.positions, sim.velocities, sim.forces  # Copies as NumPy arrays, sorted by particle ID.
plot_particles(sim, color_by="speed")
plot_history(sim)
```

### Configuration

The configuration is a dictionary with the keys of md-flexible's YAML files, or the path to such a file.
[`AllOptions.yaml`](../examples/md-flexible/input/AllOptions.yaml) lists all of them.
Both go through the same parser, so `ap.AutoPas("fallingDrop.yaml")` behaves like `md-flexible --yaml-filename
fallingDrop.yaml`, apart from what is listed under [limitations](#differences-to-md-flexible).
On top of what YAML can express, `"Objects"` can be a list of generators.

### Functors

| Class                  | Interaction | Note                                 |
|------------------------|-------------|--------------------------------------|
| `LennardJones`         | pairwise    | Auto vectorized                      |
| `LennardJonesHWY`      | pairwise    | Google Highway SIMD, see note below  |
| `LennardJonesAVX`      | pairwise    | Opt-in: `AUTOPAS_PYTHON_FUNCTOR_AVX` |
| `LennardJonesSVE`      | pairwise    | Opt-in: `AUTOPAS_PYTHON_FUNCTOR_SVE` |
| `AxilrodTellerMuto`    | triwise     | Auto vectorized                      |

> [!WARNING]
> With the SoA data layout, `LennardJonesHWY` reports a wrong potential energy and virial as soon as there are halo
> particles, as its kernels count pairs with a halo particle fully instead of half. Forces are correct.
> This is a bug in `LJFunctorHWY`, which also affects md-flexible. Use `LennardJones` if you need the energies.

Cutoff and particle properties come from the simulation.
If `functors` is omitted, the keys `functor` and `functor-3b` of the configuration decide.
A list can hold at most one native functor per interaction type, plus any number of functors written in Python.

### Writing functors in Python

```python
class Gravity(ap.PairwiseFunctor):
    def aos(self, i, j, newton3):
        dr = j.r - i.r
        force = dr / dr.dot(dr) ** 1.5
        i.add_f(force)
        if newton3:
            j.sub_f(force)

sim = ap.AutoPas(config, functors=Gravity())
```

* `aos(i, j, newton3)` has to update `i`, and `j` only if `newton3` is set. `TriwiseFunctor` has `aos(i, j, k, newton3)`.
* The particles are only valid during the call. Do not keep references to them.
* Only pairs within the cutoff reach Python. The base class takes `cutoff` to narrow this down further.
* Pass `allows_newton3=False` or `allows_non_newton3=False` to `super().__init__()` if one mode is not implemented.
* Optional overrides are `init_traversal()`, `end_traversal(newton3)`, `potential_energy()`, and `virial()`.
* With `deltaT: 0`, as in md-flexible, no halo particles are exchanged. So boundary conditions have no effect.
* An exception ends the traversal and is raised by the call to the simulation. The forces are incomplete at that point,
  so the object can not be used any further.

Python functors run with the AoS data layout only, so a simulation restricts its search space to AoS as soon as one is
involved, also for native functors of the same interaction type.
Their traversal runs on a single thread, because more threads only compete for the GIL, and each kernel call enters the
interpreter. Expect this to be orders of magnitude slower than a native functor. See [`lennardJonesCube.py`](../examples/python/lennardJonesCube.py) for a comparison.

### Generators and particle properties

`autopas.generators` holds `CubeGrid`, `CubeUniform`, `CubeGauss`, `CubeClosestPacked`, and `Sphere`.
`generate()` returns the particles as arrays, or they are handed over via `"Objects"` or `objects=`.

`ap.ParticlePropertiesLibrary(cutoff)` can replace `"Sites"` via `particle_properties=`.
`sim.particle_properties` gives access to the library of a simulation, including the mixed values.

## Differences to md-flexible

* No MPI, and no multi-site molecules. CMake rejects these configurations.
* Native functors calculate potential energy and virial by default, for the history of a simulation.
  This costs some performance. Disable it with `AUTOPAS_PYTHON_CALC_GLOBALS=OFF`.
* No progress bar, no end configuration, no output of timers and no global variable log.
* As in md-flexible, the first position update uses zero forces if none are given, see `Simulation::run()`.
  This means the energy of the very first iteration is off.

## Tests

```shell
pip install ".[test]"
pytest
```

[`environment.yml`](environment.yml) contains everything the CI uses to build and test the module.
