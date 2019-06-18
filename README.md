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
* a c++14 compiler (gcc7, clang6 and icpc 2018 are tested)

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
* to run specific tests:
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
### AutoPas and Containers
creating your container, adding particles,
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
for(auto iter = container.begin(); iter.isValid(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```
### Updating the Container
#### How
You can update the container using
```C++
ParticleContainer::updateContainer()
```
#### When it is necessary
You have to update the container when the two conditions are fullfilled:
* If you moved particles
* You want to use `iteratePairwise()` or a RegionParticleIterator

#### When it is not enough
If you moved particles by more than one interaction length.
If you are planning to move particles by a long distance,
e.g. because of boundary conditions please delete the particles and add them again:
```C++
std::vector<autopas::sph::SPHParticle> invalidParticles;
for (auto part = sphSystem.begin(); part.isValid(); ++part) {
  if (/*check*/) {
    invalidParticles.push_back(*part);
    part.deleteCurrentParticle();
  }
}
for (auto p: invalidParticles) {
  sphSystem.addParticle(p);
}
```

#### Special exceptions
* Verlet-Lists, here it is safe to not update the container as long as particles move not more than a skin radius.

## Developing AutoPas
Please look at our [contribution guidelines](https://github.com/AutoPas/AutoPas/blob/master/.github/CONTRIBUTING.md).

## Acknowledgements
* TaLPas BMBF

## Papers to cite
* TODO: Add papers
