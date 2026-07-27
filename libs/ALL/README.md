# A Load Balancing Library (ALL)

The library aims to provide an easy way to include dynamic domain-based
load balancing into particle based simulation codes. The library is
developed in the Simulation Laboratory Molecular Systems of the Jülich
Supercomputing Centre at Forschungszentrum Jülich. 

Only a brief summary is given here and more information can be found in
the [official documentation](http://slms.pages.jsc.fz-juelich.de/websites/all-website/) such as a detailed API
description, examples and further information regarding the load
balancing methods.

## Installation and Requirements

### Requirements
Base requirements:

 - C++11 capable compiler
 - MPI support
 - CMake v. 3.14 or higher

Optional requirements:

 - Fortran 2003 capable compiler (Fortran interface)
 - Fortran 2008 capable compiler (Usage of the `mpi_f08` interface)
 - VTK 7.1 or higher (Domain output)
 - Boost testing utilities
 - Doxygen and Sphinx with `breathe` (Documentation)

### Installation

 1. Clone the library from
    `https://gitlab.version.fz-juelich.de/SLMS/loadbalancing` into
    `$ALL_ROOT_DIR`.
 2. Create the build directory `$ALL_BUILD_DIR` some place else.
 3. Call `cmake -S "$ALL_ROOT_DIR" -B "$ALL_BUILD_DIR"` to set up the
    installation.  To use a specific compiler and Boost installation
    use: `CC=gcc CXX=g++ BOOST_ROOT=$BOOST_DIR cmake [...]`.
  4. To build and install the library then run: `cmake --build
     "$ALL_BUILD_DIR"` and `cmake --install "$ALL_BUILD_DIR" --prefix
     "$ALL_INSTALL_DIR"`.
     Afterwards, the built examples and library files are placed in
     `$ALL_INSTALL_DIR`.

<!-- vim: set tw=72 et ts=4 spell spelllang=en_us ft=markdown: -->
