.. _install:

Installation
============

Requirements
------------

Base requirements
 - C++11 capable compiler
 - MPI support
 - CMake v.3.14 or higher

Optional requirements
 - Fortran 2003 capable compiler (Fortran interface)
 - Fortran 2008 capable compiler (Usage of the ``mpi_f08`` interface)
 - VTK 7.1 or higher (Domain output)
 - Boost testing utilities
 - Doxygen and Sphinx with ``breathe`` (Documentation)

Installation
------------

1. Clone the library from
   ``https://gitlab.jsc.fz-juelich.de/SLMS/loadbalancing`` into
   ``$ALL_ROOT_DIR``.
2. Create the build directory ``$ALL_BUILD_DIR`` some place else.
3. Call ``cmake -S "$ALL_ROOT_DIR" -B "$ALL_BUILD_DIR"`` to set up the
   installation. Other compile time features should be set here. To use a
   specific compiler and Boost installation use: ``CC=gcc CXX=g++
   BOOST_ROOT=$BOOST_DIR cmake [...]``.
4. To build and install the library then run: ``cmake --build
   "$ALL_BUILD_DIR"`` and ``cmake --install "$ALL_BUILD_DIR"``. Afterwards,
   the built examples and library files are placed in ``$ALL_INSTALL_DIR``.

Compile time features
*********************
``-DCMAKE_INSTALL_PREFIX=$ALL_INSTALL_DIR`` (default: system dependent)
  Install into ``$ALL_INSTALL_DIR`` after compiling.
``-DCM_ALL_VTK_OUTPUT=ON`` (default: ``OFF``)
  Enable VTK output of domains. Requires VTK 7.1 or higher.
``-DCM_ALL_FORTRAN=ON`` (default: ``OFF``)
  Enable Fortran interface.  Requires Fortran 2003 capable compiler.
``-DCM_ALL_USE_F08=ON`` (default: ``OFF``)
  Enable usage of ``mpi_f08`` module for MPI. Requires Fortran 2008
  capable compiler and compatible MPI installation.
``-DCM_ALL_FORTRAN_ERROR_ABORT=ON`` (default: ``OFF``)
  Abort execution on any error when using the Fortran interface instead of
  setting the error number and leaving error handling to the user.
``-DCM_ALL_VORONOI=ON`` (default: ``OFF``)
  Enable Voronoi mesh method and subsequently compilation and linkage of
  Voro++.
``-DCM_ALL_TESTS=ON`` (default: ``OFF``)
  Turns on generation of tests. Will only generate unit tests by default
  and additional tests can be enabled with other flags, see below. Tests
  can be run from the build directory with ``ctest``. Requires the Boost
  test utilities.
``-DCM_ALL_TESTS_INTEGRATION=ON`` (default: ``OFF``)
  Enables integration/feature tests. Requires ``CM_ALL_TESTS`` to be
  enabled as well.
``-DCM_ALL_AUTO_DOC=ON`` (default: ``OFF``)
  Generate documentation using doxygen and Sphinx.
``-DCMAKE_BUILD_TYPE=Debug`` (default: ``Release``)
  Enable library internal debugging features. Using ``DebugWithOpt`` also
  turns on some optimizations.
``-DCM_ALL_DEBUG=ON`` (default: ``OFF``)
  Enable self consistency checks within the library itself.
``-DVTK_DIR=$VTK_DIR`` (default: system dependent)
  Path to an explicit VTK installation. Make sure that VTK is compiled
  with the same MPI as the library and subsequently your code.

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
