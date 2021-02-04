# ![AutoPas](https://raw.githubusercontent.com/AutoPas/AutoPas/master/docs/graphics/AutoPasLogo_Large.svg "Title")

AutoPas is a node-level auto-tuned particle simulation library developed
in the context of the [**TaLPas**](http://www.talpas.de) project.
[![Build Status](https://www5.in.tum.de/jenkins/mardyn/buildStatus/icon?job=AutoPas-Multibranch/master)](https://www5.in.tum.de/jenkins/mardyn/job/AutoPas-Multibranch/job/master/)

## Documentation
The documentation can be found at our website:
 <https://autopas.github.io/doxygen_documentation/git-master/>

Alternatively you can build the documentation on your own:
* requirements: [Doxygen](http://www.doxygen.nl/)
* `make doc_doxygen`

## Requirements
* CMake 3.14 or newer
* make (build-essentials) or ninja
* a c++17 compiler (gcc7, clang8 and ~~icpc 2019~~ are tested.)

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

AutoPas relies on a small number of dependencies. By default AutoPas looks for
installed versions of those libraries but it can also be forced to (selectively)
use bundled versions. To make use of this feature, call `cmake` with:
```bash
cmake -D spdlog_ForceBundled=ON    # replace spdlog by the lib you want to force
```
Or better have a look at the variables exposed in `ccmake`. 

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
Many IDEs (e.g., CLion) have integrated support for googletest and you can debug the tests directly within the IDE.

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
For that we provide some basic Particle classes defined in `src/autopas/molecularDynamics` or `src/autopas/sph` 
that you can use either directly or you can write your own Particle class by inheriting from one of the provided
classes or from `autopas::Particle`.

Important parts to implement:
* `enum AttributeNames`
* Definition of a matching `SoAArraysType`
* Getter and setter connecting the `AttributeNames` and actual members.

### Custom Functors
Once you have defined your particle you can start with the functor class.

#### Definition
Importatnt parts to implement:
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
* owned: Particles that belong to this AutoPas instance. These particles are either inside of the boundary of the AutoPas instance or very close to the boundary (less than a distance of skin/2 away). If a particle is added via `addParticle()`, it is automatically added as an owned particle. An owned particle can explicitly be removed by deleting the particle using an iterator (`autoPas.deleteParticle(iterator)`). On an update of the AutoPas container (using `updateContainer()`) owned particles that move outside of the boundary of its parent AutoPas container are returned.
* halo: Particles that do not belong to the current AutoPas instance. These normally are ghost particles arising from either periodic boundary conditions or particles of a neighboring AutoPas object (if you split the entire domain over multiple AutoPas objects, i.e., you use a domain decomposition algorithm). The halo particles are needed for the correct calculation of the pairwise forces. On update of the AutoPas container, halo particles are deleted (note that not every call to `updateContainer()` does this!, see [Simulation Loop](#simulation-loop)).
* dummy: Particles that are about to be deleted or that act as filler for certain algorithms. These particles do not affect the force calculation.

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
The default parameter is `haloAndOwned`, which is also used for range-based for loops.

Analogously to `begin()`, `cbegin()` is also defined, which guarantees to return a `const_iterator`.

### Simulation Loop
One simulation loop should always consist of the following phases:

1. Updating the Container:
   ```cpp
   auto [invalidParticles, updated] = autoPas.updateContainer();
   ```
   This call will potentially trigger an update of the container inside of AutoPas. The update will be performed if either
   
    a. The AutoTuner collected enough samples for the current configuration and will move to the next one OR
    
    b. The rebuild frequency of the container is reached.
   
   If the update is performed, the returned bool `updated` is `true`. The returned vector `invalidParticles` consists of the particles that are not anymore within the boundaries of this container and hence are deleted from it. These are particles that were previously owned by this AutoPas container but have left the boundary of this container, i.e., their current position resides outside of the container.
   
   If the update is not performed, `updated` will be false and the returned vector `invalidParticles` will be empty.
   An update is sometimes skipped to ensure that containers do not change, which allows containers to reuse neighbor lists thus enabling better performance.

2. Handling the leaving particles
   * This step can be skipped if `updated` was false. If you use multiple MPI instances, you have to ensure that all instances rebuild during the same time step. This is guaranteed if the sampling frequency is the same as (or a multiple of) the rebuild frequency.
   * Apply boundary conditions on them
   * Potentially send them to other mpi-processes, skip this if MPI is not needed
   * Add them to the containers using
      ```cpp
      autoPas.addParticle(particle)
      ```

3. Handle halo particles:
   * This step always has to be performed, even if `updated` was false.
   * Identify the halo particles by use of AutoPas' iterators and send them in a similar way as the leaving particles.
   * Add the particles as haloParticles using
      ```cpp
      autoPas.addOrUpdateHaloParticle(haloParticle)
      ```

4. Perform an iteratePairwise step.
   ```cpp
   autoPas.iteratePairwise(functor);
   ```

### Inserting additional particles
Before inserting additional particles (e.g. through a grand-canonical thermostat ),
you always have to enforce a containerUpdate on ALL AutoPas instances, i.e.,
on all mpi processes, by calling
```cpp
autoPas.updateContainerForced();
```
This will invalidate the internal neighbor lists and containers.

### Using multiple functors

AutoPas is able to work with simulation setups using multiple functors that describe different forces.
A good demonstration for that is the sph example found under examples/sph or examples/sph-mpi.
There exist some things you have to be careful about when using multiple functors:
* If you use multiple functors it is necessary that all functors support the same newton3 options. If there is one functor not supporting newton3, you have to disable newton3 support for AutoPas by calling
  ```cpp
  autoPas.setAllowedNewton3Options({false});
  ```
* If you have `n` functors within one iteration and update the particle position only at the end or start of the iteration, the rebuildFrequency and the samplingRate have to be a multiple of `n`.
* Functors must be marked as (not) relevant for tuning by specifying `Functor::isRelevantForTuning()`. Functors marked as relevant should have a near-identical performance profile otherwise the sampling of configurations will be distorted. It is recommended, to only mark the most expensive functor as relevant.

## Developing AutoPas
Please look at our [contribution guidelines](https://github.com/AutoPas/AutoPas/blob/master/.github/CONTRIBUTING.md).

## Acknowledgements
This work was financially supported by:
* the Federal Ministry of Education and Research, Germany, project “Task-based load balancing and auto-tuning in particle simulations” (TaLPas) 8 , grant numbers 01IH16008A and 01IH16008B.

## Papers to cite
* F. A. Gratl, S. Seckler, N. Tchipev, H.-J. Bungartz and P. Neumann: [AutoPas: Auto-Tuning for Particle Simulations](https://ieeexplore.ieee.org/document/8778280) [BibTeX](https://mediatum.ub.tum.de/export/1535848/bibtex) [MediaTUM](https://mediatum.ub.tum.de/1535848), In 2019 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), Rio de Janeiro, May 2019.
* S. Seckler, F. Gratl, M. Heinen, J. Vrabec, H.-J. Bungartz, P. Neumann: [AutoPas in ls1 mardyn: Massively parallel particle simulations with node-level auto-tuning](https://www.sciencedirect.com/science/article/abs/pii/S1877750320305901) [BibTeX](https://mediatum.ub.tum.de/export/1595680/bibtex) [MediaTUM](https://mediatum.ub.tum.de/1595680), In Journal of Computational Science, Volume 50, 2021.
