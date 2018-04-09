# Autopas
Autopas is a node-level auto-tuned particle simulation library developed
in the context of the **TaLPas** project.

## Documentation
The documentation can be found at our website:
 https://www5.in.tum.de/autopas/doxygen_doc/html/

Alternatively you can build the documentation on your own:
* requirements:
 doxygen
* `make doc_doxygen`


## Requirements
* cmake 3.3 or newer
* build-essentials (make)
* a c++11 compiler


## Building autopas
build instructions:
```
mkdir build
cd build
cmake ..
make
```


## Testing
to run tests:
```
make test
```
or using the ctest environment:
```
ctest
```
to get verbose output:
```
ctest --verbose
```
* to run specific tests:
use the --gtest_filter variable:
```
./tests/testAutopas/runTests --gtest_filter=ArrayMathTest.testAdd*
```
or use the GTEST_FILTER environment variable:
```
GTEST_FILTER="ArrayMathTest.testAdd*" ctest --verbose
```


## Examples
As AutoPas is only library for particle simulations it itself is not able to run simulations.
We have, however, included a variety of examples in the **examples** directory. The examples include:
* Molecular dynamics simulations with 1 centered Lennard-Jones particles.
* Smoothed particle hydrodynamics simulations
* Gravity simulations


## Using Autopas

Steps to using Autopas in your particle simulation program:
1. **Defining a Custom Particle Class** <br/>
First you will need to define a particle class.
For that we provide some basic Particle classes defined
in `src/particles/` that you can use either directly
or you can write your own Particle class by inheriting from
one of the provided classes.
```
class SPHParticle : public autopas::Particle {

}
```
2. **AutoPas and Containers**<br>
adding particles, iterating through containers
3. **Defining Functors** <br/>
Once you have defined your particle you can start
4. **Iterating Through Particles**



## Developing Autopas
* We use google code style.
* code style can be build with `make clangformat`
* requirements:
	clang-format

## Acknowledgements
* TaLPas BMBF

## Papers to cite
* TODO: Add papers
