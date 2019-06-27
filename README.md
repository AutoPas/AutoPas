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
* cmake 3.13 or newer
* make (build-essentials) or ninja
* a c++17 compiler (gcc7, clang8 and icpc 2018 are tested)

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
#### Usage

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
### Simulation Loop
One simulation loop should always consist of the following phases:

1. Updating the Container, which returns a vector of all invalid == leaving particles!
   ```C++
   auto invalidParticles = autoPas.updateContainer();
   ```

2. Handling the leaving particles
   * Apply boundary conditions on them
   
   * Potentially send them to other mpi-processes, skip this if MPI is not needed
   
   * Add them to the containers using
      ```C++
      autoPas.addParticle(particle)
      ```

3. Handle halo particles:
   * Identify the halo particles by use of AutoPas' iterators and send them in a similar way as the leaving particles.

   * Add the particles as haloParticles using 
      ```C++
      autoPas.addOrUpdateHaloParticle(haloParticle)
      ```

4. Perform an iteratePairwise step.
   ```C++
   autoPas.iteratePairwise(functor);
   ```

In some iterations step 1. will return an empty list of invalid particles to benefit of not rebuilding the containers and the associated neighbor lists.

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
* TaLPas BMBF

## Papers to cite
* TODO: Add papers
