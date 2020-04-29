# AutoPas
AutoPas is a node-level auto-tuned particle simulation library developed
in the context of the **TaLPas** project. [![Build Status](https://www5.in.tum.de/jenkins/mardyn/buildStatus/icon?job=AutoPas-Multibranch/master)](https://www5.in.tum.de/jenkins/mardyn/job/AutoPas-Multibranch/job/master/)

## Documentation
The documentation can be found at our website:
 <https://www5.in.tum.de/AutoPas/doxygen_doc/master/>

Alternatively you can build the documentation on your own:
* requirements: [Doxygen](http://www.doxygen.nl/)
* `make doc_doxygen`

## Requirements
* cmake 3.14 or newer
* make (build-essentials) or ninja
* a c++17 compiler (gcc7, clang8 and icpc 2019 are tested)

## Building AutoPas
build instructions for make:
```bash
mkdir build
cd build
cmake ..
make
```
if you want to use another compiler, specify it at the first cmake call, e.g.:
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

## Testing
### Running Tests
to run tests:
```bash
make test
# or
ninja test
```
or using the ctest environment:
```bash
ctest
```
to get verbose output:
```bash
ctest --verbose
```
#### How to run specific tests

use the --gtest_filter variable:
```bash
./tests/testAutopas/runTests --gtest_filter=ArrayMathTest.testAdd*
```
or use the GTEST_FILTER environment variable:
```bash
GTEST_FILTER="ArrayMathTest.testAdd*" ctest --verbose
```
or `ctest` arguments like `-R` (run tests matching regex) and `-E` (exclude tests matching regex)
```bash
ctest -R 'Array.*testAdd' -E 'Double'
```

### Debugging Tests
Find out the command to start your desired test with `-N` aka. `--show-only`:
```bash
ctest -R 'Array.*testAdd' -N
```
Start the test with `gdb`
```bash
gdb --args ${TestCommand}
```

## Examples
As AutoPas is only a library for particle simulations it itself is not able to run simulations.
We have, however, included a variety of examples in the **examples** directory. The examples include:
* Molecular dynamics simulations with 1 centered Lennard-Jones particles.
* Smoothed particle hydrodynamics simulations
* Gravity simulations

## Using AutoPas
Steps to using AutoPas in your particle simulation program:

### Defining a Custom Particle Class
First you will need to define a particle class.
For that we provide some basic Particle classes defined
in `src/particles/` that you can use either directly
or you can write your own Particle class by inheriting from
one of the provided classes.
```C++
class SPHParticle : public autopas::Particle {

}
```

### Functors
Once you have defined your particle you can start with functors;
#### Definition
TODO
#### Usage
TODO

### Particle Ownership
An AutoPas container normally saves two different types of particles:
* owned particles: Particles that belong to the AutoPas instance. These particles lie either inside of the boundary of the AutoPas instance or very close to the boundary (less than a distance of skin/2 away). If a particle is added via `addParticle()`, it is automatically added as owned particle. An owned particle can explicitly be removed by deleting the particle using an iterator (`autoPas.deleteParticle(iterator)`). On an update of the AutoPas container (using `updateContainer()`) owned particles that move outside of the boundary of its parent AutoPas container are returned.
* halo particles: Particles that do not belong to the current AutoPas instance. These normally are ghost particles arising from either periodic boundary conditions or particles of a neighboring AutoPas object (if you split the entire domain over multiple AutoPas objects, i.e., you use a domain decompositioning algorithm). The halo particles are needed for the correct calculation of the pairwise forces. On an update of the AutoPas container, halo particles are deleted (note that not every call to `updateContainer()` does this!, see [Simulation Loop]).

### Iterating Through Particles
Iterators to iterate over particle are provided.
The particle can be accesses using `iter->` (`*iter` is also possible).
When created inside a OpenMP parallel region, work is automatically spread over all iterators.
```C++
#pragma omp parallel
for(auto iter = autoPas.begin(); iter.isValid(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```
For convenience the `end()` method is also implemented for the AutoPas class.
```C++
#pragma omp parallel
for(auto iter = autoPas.begin(); iter != autoPas.end(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```
You might also use range-based for loops:
```C++
#pragma omp parallel
for(auto& particle : autoPas) {
  // user code:
  auto position = particle.getR();
}
```

To iterate over a subset of particles, the `getRegionIterator(lowCorner, highCorner)`
method can be used:
```C++
#pragma omp parallel
for(auto iter = autoPas.getRegionIterator(lowCorner, highCorner); iter != autoPas.end(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```

Both `begin()` and `getRegionIterator()` can also take the additional parameter `IteratorBehavior`,
which indicates over which particles the iteration should be performed.
```C++
enum IteratorBehavior {
  haloOnly,     /// iterate only over halo particles
  ownedOnly,    /// iterate only over owned particles
  haloAndOwned  /// iterate over both halo and owned particles
};
```
The default parameter is `haloAndOwned`, which is also used for range-based for loops.

Analogously to `begin()`, `cbegin()` is also defined, which guarantees to return a
`const_iterator`.

### Simulation Loop
One simulation loop should always consist of the following phases:

1. Updating the Container:
   ```C++
   auto [invalidParticles, updated] = autoPas.updateContainer();
   ```
   This call will potentially trigger an update of the container inside of AutoPas. The update will be performed if either
   
    a. Enough samples are collected for the current configuration OR
    
    b. The rebuild frequeny is reached.
   
   If the update is performed, the returned bool `updated` is `true` and the returned vector `invalidParticles` consists of the particles that are leaving this container. These are particles which were previously owned by this AutoPas container, but have left the boundary of this container, i.e., their current position resides outside of the container.
   
   If the update is not performed, `updated` will be false and the returned vector `invalidParticles` will be empty.
   An update is sometimes skipped to ensure that containers do not change and to enable the best possible performance for containers that use neighbor lists.

2. Handling the leaving particles
   * This step can be skipped if `updated` was false. If you use multiple MPI instances, you have to ensure that all instances rebuild at the same time steps. This is guaranteed, if the sampling frequency is the same as (or a multiple of) the rebuild frequency.
   
   * Apply boundary conditions on them

   * Potentially send them to other mpi-processes, skip this if MPI is not needed

   * Add them to the containers using
      ```C++
      autoPas.addParticle(particle)
      ```

3. Handle halo particles:
   * This step always has to be performed, even if `updated` was false.
   
   * Identify the halo particles by use of AutoPas' iterators and send them in a similar way as the leaving particles.

   * Add the particles as haloParticles using
      ```C++
      autoPas.addOrUpdateHaloParticle(haloParticle)
      ```

4. Perform an iteratePairwise step.
   ```C++
   autoPas.iteratePairwise(functor);
   ```


### Using multiple functors

AutoPas is able to work with simulation setups using multiple functors that describe different forces.
A good demonstration for that is the sph example found under examples/sph or examples/sph-mpi.
There exist some things you have to be careful about when using multiple functors:
* If you use multiple functors it is necessary that all functors support the same newton3 options. If there is one functor not supporting newton3, you have to disable newton3 support for AutoPas by calling
  ```C++
  autoPas.setAllowedNewton3Options({false});
  ```

* If you have `n` functors within one iteration and update the particle position only at the end or start of the iteration, the rebuildFrequency and the samplingRate have to be a multiple of `n`.

### Inserting additional particles
Before inserting additional particles (e.g. through a grand-canonical thermostat ),
you always have to enforce a containerUpdate on ALL AutoPas instances, i.e.,
on all mpi processes, by calling
```C++
autoPas.updateContainerForced();
```
This will invalidate the internal neighbor lists and containers.

## Developing AutoPas
Please look at our [contribution guidelines](https://github.com/AutoPas/AutoPas/blob/master/.github/CONTRIBUTING.md).

## Acknowledgements
This work was financially supported by:
* the Federal Ministry of Education and Research, Germany, project “Task-based load balancing and auto-tuning in particle simulations” (TaLPas) 8 , grant numbers 01IH16008A and 01IH16008B.

## Papers to cite
* F. A. Gratl, S. Seckler, N. Tchipev, H.-J. Bungartz and P. Neumann: [AutoPas: Auto-Tuning for Particle Simulations](https://ieeexplore.ieee.org/document/8778280) [BibTeX](https://mediatum.ub.tum.de/services/export/node/1535848/?format=template_test&mask=bibtex&lang=de&template=$$[defaultexport]$$&mimetype=text/plain) [MediaTum](https://mediatum.ub.tum.de/1535848), In 2019 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), Rio de Janeiro, May 2019.

