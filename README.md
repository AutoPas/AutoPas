
# ![AutoPas](https://raw.githubusercontent.com/AutoPas/AutoPas/master/docs/graphics/AutoPasLogo_Large.svg "Title")

AutoPas is a node-level auto-tuned particle simulation library developed
in the context of the [**TaLPas**](http://www.talpas.de) project.
[![CI Status](https://github.com/AutoPas/AutoPas/actions/workflows/TestSuites.yaml/badge.svg)](https://github.com/AutoPas/AutoPas/actions/workflows/TestSuites.yaml)

## Documentation
The documentation can be found at our website:
 <https://autopas.github.io/doxygen_documentation/git-master/>

Alternatively you can build the documentation on your own:
* requirements: [Doxygen](http://www.doxygen.nl/)
* `make doc_doxygen`

## Requirements
* CMake 3.14 or newer
* make (build-essentials) or ninja
* a C++17 compiler (gcc11, clang13, and ~~icpc 2019~~ are tested.)

## Building AutoPas
build instructions for make:
```bash
mkdir build
cd build
cmake ..
make
```
if you want to use a specific compiler, specify it at the first CMake call, e.g.:
```bash
mkdir build
cd build
CC=clang CXX=clang++ cmake ..
make
```
if you would like to use ninja instead of make:
```bash
mkdir build
cd build
cmake -G Ninja ..
ninja
```
### Building AutoPas on a Cluster
HPC clusters often use module systems. CMake is sometimes not able to
correctly detect the compiler you wished to use. If a wrong compiler is
found please specify the compiler explicitly, e.g. for gcc:
```bash
mkdir build
cd build
CC=`which gcc` CXX=`which g++` cmake ..
make
```

AutoPas relies on a small number of dependencies. By default, AutoPas looks for
installed versions of those libraries, but it can also be forced to (selectively)
use bundled versions. To make use of this feature, call `cmake` with:
```bash
cmake -D spdlog_ForceBundled=ON    # replace spdlog by the lib you want to force
```
Or better, have a look at the variables exposed in `ccmake`. 

## Testing

AutoPas uses [googletest](https://github.com/google/googletest) as testing 
framework and exposes tests to ctest, the CMake test driver.

### Running Tests

There are multiple possibilities. In order of recommendation:

1. Using `ctest`:
   ```bash
   ctest # add --verbose for more details on the tests
   ```
   To only run specific tests use arguments like `-R` (run tests matching regex) and `-E` (exclude tests matching regex)
   ```bash
   ctest -R 'Array.*testAdd' -E 'Double'
   ```
2. Using the `make` target:
   ```bash
   make test
   ```
3. Directly launching the test executable:
   ```bash
   tests/testAutopas/runTests
   ```
   To only run specific tests use arguments
   ```bash
   tests/testAutopas/runTests --gtest_filter=ArrayMathTest.testAdd*
   ```
 
### Debugging Tests
Many IDEs (e.g., CLion) have integrated support for googletest, and you can debug the tests directly within the IDE.

If you prefer `gdb`:
1. Find out the command to start your desired test with `-N` aka. `--show-only`:
   ```bash
   ctest -R 'Array.*testAdd' -N
   ```
2. Start the test with `gdb`
   ```bash
   gdb --args ${TestCommand}
   ```

## Examples
As AutoPas is only a library, it is not able to run simulations by itself.
We have, however, included a few example proxy applications in the **examples** directory.
The examples include:
* [md-flexible](examples/md-flexible): Molecular dynamics simulations with single centered Lennard-Jones particles.
* Smoothed particle hydrodynamics simulations

## Using AutoPas
Steps to using AutoPas in your particle simulation program:

### Custom Particles
First you will need to define a particle class which will be passed to AutoPas as template Argument.
For that we provide some basic Particle classes defined in `applicationLibrary/molecularDynamicsLibrary` or `src/autopas/sph` 
that you can use either directly or you can write your own Particle class by inheriting from one of the provided
classes or from `autopas::Particle`.

Important parts to implement:
* `enum AttributeNames`
* Definition of a matching `SoAArraysType`
* Getter and setter connecting the `AttributeNames` and actual members.

### Custom Functors
Once you have defined your particle you can start with the functor class.

#### Definition
Important parts to implement:
* Actual force calculations: `AoSFunctor()` and all Versions of `SoAFunctor*()` 
* Newton3 characteristics of the force: `allowsNewton3()` and `allowsNonNewton3()`
* Input and output variables of the force calculation via: `getComputedAttr()` and `getNeededAttr()`

#### Usage
Each functor is applied to AutoPas via:
```bash
autoPas.iteratePairwise(&myFunctor);
```

### Particle Ownership
Particles saved in an AutoPas container can be one of two possible states:
* owned: Particles that belong to this AutoPas instance. 
  These particles are typically inside the boundary of the AutoPas instance.
  If a particle is added via `addParticle()`, it is automatically added as an owned particle.
  An owned particle can explicitly be removed by deleting the particle using an iterator (`autoPas.deleteParticle(iterator)`).
  On an update of the AutoPas container (using `updateContainer()`) owned particles that moved outside the boundary of its parent AutoPas container are returned.
* halo: Particles that do not belong to the current AutoPas instance.
  These normally are ghost particles arising from either periodic boundary conditions or particles of a neighboring AutoPas object
  (if you split the entire domain over multiple AutoPas objects, i.e., you use a domain decomposition algorithm).
  The halo particles are needed for the correct calculation of the pairwise forces.
  On update of the AutoPas container, halo particles are deleted (see <a href="#simulation-loop">Simulation Loop</a>).
* dummy: Particles that are deleted or that act as filler for certain algorithms. These particles do not affect the force calculation.

### Iterating Through Particles
Iterators to iterate over particle are provided.
The particle can be accesses using `iter->` or `*iter`.
When created inside a OpenMP parallel region, work is automatically spread over all threads.
```cpp
#pragma omp parallel
for(auto iter = autoPas.begin(); iter.isValid(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```
For convenience the `end()` method is also implemented for the AutoPas class so you might also use range-based for loops:
```cpp
#pragma omp parallel
for(auto& particle : autoPas) {
  // user code:
  auto position = particle.getR();
}
```

To iterate over a subset of particles, the `getRegionIterator(lowCorner, highCorner)` method can be used:
```cpp
#pragma omp parallel
for(auto iter = autoPas.getRegionIterator(lowCorner, highCorner); iter != autoPas.end(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```

Both `begin()` and `getRegionIterator()` can also take the additional parameter `IteratorBehavior`,
which indicates over which particles the iteration should be performed. See [autopas::IteratorBehavior
](https://autopas.github.io/doxygen_documentation/git-master/namespaceautopas.html#a520fefd51e4555074cd16e7c3fd19c42) for possible options and details.
The default parameter is `ownedOrHalo`, which is also used for range-based for loops.

Analogously to `begin()`, `cbegin()` is also defined, which guarantees to return a `const_iterator`.

Iterators are not guaranteed to be valid after particle insertion. 
However, particle deletion while iterating is supported via `autoPas.deleteParticle(iterator)`. 
After deletion the `++` operator has to be called:
```cpp
#pragma omp parallel
for(auto iter = autoPas.getIterator(); iter != autoPas.end(); ++iter) {
  autoPas.deleteParticle(iterator);
}
```

### Logging

AutoPas contains multiple loggers with different purposes that can help to shed light into the black box.
Under the hood, they use [spdlog](https://github.com/gabime/spdlog). 
When deactivated via `CMake` these loggers do not add any run time overhead. 

#### AutoPasLog
This is the main, general purpose logger. It supports all spdlog-levels. These levels can be (de)activated at compile
time via the `CMake` variable `AUTOPAS_MIN_LOG_LVL`. At run time, this logger's compiled levels can be set e.g. via: 
`autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);`

At debug level, this logger will print the full configuration of every call to `autopas::AutoPas::iteratePairwise()`.

#### GaussianClusterLogger
Creates a graph representation of the Gaussian cluster model that was created during the simulation.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_GAUSSIANCLUSTER`.

#### IterationLogger
Creates a csv file containing information about the configuration and timings of every single pairwise iteration.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_ITERATIONS`.

#### OctreeLogger
Creates a vtk file to visualize the Octree particle container.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_OCTREE`.

#### PredictionLogger
Creates a csv containing the predictions made by the PredictiveTuning strategy.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_PREDICTIONS`.

#### TuningDataLogger
Creates a csv containing all data that is collected for tuning purposes. This is the raw data that is available to 
the tuning algorithms. 
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_TUNINGDATA`.

#### TuningResultLogger
Creates a csv containing the results of every tuning phase. Useful if only the high level end results are of interest.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_TUNINGRESULTS`.

### Simulation Loop
*TODO SHOW WHOLE LOOP WITH EXAMPLE!*
One simulation loop should always consist of the following phases:

1. Updating the Container:
   ```cpp
   auto invalidParticles = autoPas.updateContainer();
   ```
   This call will trigger an update of the container inside AutoPas. 
   The returned vector `invalidParticles` consists of the particles that were previously owned by this AutoPas container
   but have left the boundary of this container, i.e., their current position resides outside the container.

2. Handling the leaving particles
   * Apply boundary conditions on them
   * Potentially send them to other mpi-processes, skip this if MPI is not needed
   * Add them to the containers using
      ```cpp
      autoPas.addParticle(particle)
      ```

3. Handle halo particles:
   * Identify the halo particles by use of AutoPas' iterators and send them in a similar way as the leaving particles.
   * Add the particles as haloParticles using
      ```cpp
      autoPas.addHaloParticle(haloParticle)
      ```

4. Perform an iteratePairwise step.
   ```cpp
   autoPas.iteratePairwise(functor);
   ```

### Inserting additional particles
Additional particles (e.g. through a grand-canonical thermostat), can be inserted at any point in the simulation loop.
For periodic boundary conditions, or in an MPI-parallel simulation, you, as the user, is responsible for inserting the appropriate halo particles.

### Internal Verlet-like container behavior
The behavior described in this section is normally opaque to users of AutoPas. The only exception to this rule is that
particles should not be moved more than skin/2 within the specified Verlet rebuild frequency. This restriction is due to
the internally used Verlet-like container behavior in which the actual container is not updated in every time step and
particles are not necessarily sorted into the correct cells. This allows the reuse of neighbor lists throughout multiple
time steps and is necessary for a performant implementation of our Verlet containers.

We do, however, still provide a linked cells-like interface to a user of AutoPas, i.e., a container appears to be
updated every time step, leaving particles are returned at every time step and particles can be deleted and added
independently to the internal state of the container. Internally we make this possible, by using partial container
updates, which collect leaving particles while marking left particles and halo particles as dummy. Additionally, we
maintain a particle buffer that allows to add particles to AutoPas without modifying the underlying container. This
particle buffer is considered in the force calculation and when iterating through particles.

Another performance optimization is made possible by allowing to reuse the neighbor list entries of halo particles of
previous time steps. While the actual particles have already been implicitly deleted (marked as dummy), they still
exist. For their reuse, we try to add halo particles in their original memory location. If that is, however, not
possible, we add them to another particle buffer (the haloParticleBuffer).

Additional information can be found in [PR 642](https://github.com/AutoPas/AutoPas/pull/642)

### Using multiple functors

AutoPas is able to work with simulation setups using multiple functors that describe different forces.
A good demonstration for that is the sph example found under examples/sph or examples/sph-mpi.
There exist some things you have to be careful about when using multiple functors:
* If you use multiple functors it is necessary that all functors support the same newton3 options.
  If there is one functor not supporting newton3, you have to disable newton3 support for AutoPas by calling
  ```cpp
  autoPas.setAllowedNewton3Options({false});
  ```
* If you have `n` functors within one iteration and update the particle position only at the end or start of the iteration,
  the rebuildFrequency and the samplingRate have to be a multiple of `n`.
* Functors must be marked as (not) relevant for tuning by specifying `Functor::isRelevantForTuning()`.
  Functors marked as relevant should have a near-identical performance profile otherwise the sampling of configurations will be distorted.
  It is recommended, to only mark the most expensive functor as relevant.

## Developing AutoPas
Please look at our [contribution guidelines](https://github.com/AutoPas/AutoPas/blob/master/.github/CONTRIBUTING.md).

For profiling the compile-time, the `cmake` option `AUTOPAS_COMPILE_TIME_PROFILING` can be turned on. This enables gcc's -`ftime-report` and clang's `-ftime-trace`. 
It is recommended to use clang, as its output is more detailed.
`-ftime-trace` generates a .json file for each compilation unit next to the generated object file (inside one of the CMakeFiles directories).
Chrome has a built-in tool for viewing these files in a flame graph. It can be accessed through the URL `chrome://tracing`.

## Acknowledgements
This work was financially supported by:
* the Federal Ministry of Education and Research, Germany, project “Task-based load balancing and auto-tuning in particle simulations” (TaLPas) 8 , grant numbers 01IH16008A and 01IH16008B.

## Papers to cite
* F. A. Gratl, S. Seckler, H.-J. Bungartz and P. Neumann: [N Ways to Simulate Short-Range Particle Systems: Automated Algorithm Selection with the Node-Level Library AutoPas](https://www.sciencedirect.com/science/article/abs/pii/S001046552100374X), In Computer Physics Communications, Volume 273, 2022. ([BibTeX](https://mediatum.ub.tum.de/export/1638766/bibtex), [MediaTUM](https://mediatum.ub.tum.de/1638766))
* F. A. Gratl, S. Seckler, N. Tchipev, H.-J. Bungartz and P. Neumann: [AutoPas: Auto-Tuning for Particle Simulations](https://ieeexplore.ieee.org/document/8778280), In 2019 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), Rio de Janeiro, May 2019. ([BibTeX](https://mediatum.ub.tum.de/export/1535848/bibtex), [MediaTUM](https://mediatum.ub.tum.de/1535848))
* S. Seckler, F. Gratl, M. Heinen, J. Vrabec, H.-J. Bungartz, P. Neumann: [AutoPas in ls1 mardyn: Massively parallel particle simulations with node-level auto-tuning](https://www.sciencedirect.com/science/article/abs/pii/S1877750320305901), In Journal of Computational Science, Volume 50, 2021. ([BibTeX](https://mediatum.ub.tum.de/export/1595680/bibtex), [MediaTUM](https://mediatum.ub.tum.de/1595680))
