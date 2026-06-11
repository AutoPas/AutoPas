.. _examples:

Examples
========
In the ``/examples`` sub directory several C++ and Fortran example
programs are place, as well as example projects. In some cases, where
`ALL` would be confusing, the library is referred to as `libALL`.

Projects
--------
To show how to integrate our library into your project, there are example
CMake and GNU make projects.

CMake
^^^^^
In the ``/examples`` directory are two example CMake projects. One uses a
separately installed (or at least compiled) ALL, and the other includes
the library as a subdirectory and builds it along with the main project.

When copying the example projects, make sure the symbolic links are still
resolved. The example programs of the project are just linked from the
examples, with the exception of ``ALL_test.cpp``. This program is
modified, so different feature flags are used for VTK and Voronoi, to test
external usage of these features. The change happens in the shell script
building the project.

In both cases, if VTK output is enabled, CMake must be able to find it. If
it is not able to by default, or a specific version should be used, set
VTK_DIR to the directory containing VTK's CMake configuration. For our
test system, which uses VTK 7.1, this is the ``/lib/cmake/vtk-7.1``
subdirectory of the install prefix of VTK.


Package
"""""""
The example project using ``find_package`` is available in
``/example/CMakeProject``. The full build steps are visible from
``build_all.sh``, which first compiles the library and then the example
project, after setting ``ALL_DIR`` (and ``VTK_DIR``).

In general, it suffices to just use ``find_package(ALL 0.9)`` in your
project's ``CMakeLists.txt`` and then link the corresponding executables
against this using ``target_link_library(YOUR_TARGET PUBLIC ALL:ALL)``.
The libraries targets are exported in the ``ALL::`` namespace. The Fortran
module is then ``ALL::ALL_fortran``.


Subdirectory
""""""""""""
This project includes the library directly in the source tree via
``add_subdirectory``. To avoid cyclic symbolic links, this version does
not run out of the box, since we need the source tree of the library as a
subdirectory of the project's directory. So copy the contents of
``/example/CMakeProjectSubdir``. The files assume the library to reside in
the subdirectory ``all`` and, if VTK is enabled, the VTK installation
in ``vtk_bin``. Also remember, that the symbolic links still resolve. Then
you can just run CMake and the project, along with the library, will be
build. This is also referenced in ``build_all.sh``.

The targets you have to link your executables against are, however, ``ALL``
or ``ALL_fortran`` respectively. The library is only namespaced if using
the aforementioned ``find_package`` method. And some name collisions may
therefore also occur. Using this method causes the library to be
automatically be build, but also to be rebuild every time the build
directory is cleaned.


GNU Make
^^^^^^^^
An example of integrating ALL into a traditional GNU Makefile project is
shown in ``/example/MakefileProject``. A symbolic link to the VTK
installation directory is assumed to exist as ``vtk_bin``, otherwise the
Makefile must be updated (see ``VTK_DIR``). It also assumes its directory
in the source tree as relative position to the ALL source directory.
Moving the project needs updating of ``ALL_SOURCE``. It will automatically
compile ALL using the options provided in ``LIBALL_CONFIGURE``.


Programs
--------

``ALL_test``
^^^^^^^^^^^^
MPI C++ code that generates a particle distribution over the domains. At
program start the domains form an orthogonal decomposition of the cubic 3d
system. Each domain has the same volume as each other domain. Depending
on the cartesian coordinates of the domain in this decomposition, a number
of points is created on each domain. The points are then distributed
uniformly over the domain. For the number of points a polynomial formula
is used to create a non-uniform point distribution. As an estimation of
the work load in this program the number of points within a domain was
chosen. There is no particle movement included in the current version of
the example code, therefore particle communication only takes place due to
the shift of domain borders.

The program creates three output files in its basic version:

``minmax.dat``
  Columns
    ``<iteration count> <W_min/W_avg> <W_max_W_avg>
    <(W_max-W_min)/(W_max+W_min)>``
  Explanation
    In order to try to give a comparable metric for the success of the
    load balancing procedure the relative difference of the minimum and
    maximum work loads in the system in relation to the average work load
    of the system are given. The last column gives an indicator for the
    inequality of work in the system, in a perfectly balanced system, this
    value should be equal to zero.
    
``stddev.txt``
  Columns
    ``<iteration count> <std_dev>``
  Explanation
    The standard deviation of work load over all domains.

``domain_data.dat``
  Columns
    ``<rank> <x_min> <x_max> <y_min> <y_max> <z_min> <z_max> <W>``
  Explanation
    There are two blocks of rows, looking like this, the first block is
    the starting configuration of the domains, which should be a uniform
    grid of orthogonal domains. The second block is the configuration
    after 200 load balancing steps. In each line the MPI rank of the
    domain and the extension in each cartesian direction is given, as well
    as the work load (the number of points).

If the library is compiled with VTK support, ``ALL_test`` also creates a
VTK based output of the domains in each iteration step. This output will
be placed in a ``vtk_outline`` sub directory, which is created if it does
not already exist. The resulting output can be visualized with tools like
ParaView or VisIt.

``ALL_test_f`` and ``ALL_test_f_obj``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Fortran example provides a more basic version of the test program
``ALL_test``, as its main goal is to show the functionality of the Fortran
interface. The code creates a basic uniform orthogonal domain
decomposition and creates an inhomogeneous particle distribution over
these. Only one load balancing step is executed and the program prints
out the domain distribution of the start configuration and of the final
configuration.

``ALL_Staggered`` and ``ALL_Staggered_f``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These create a very simple load balancing situation with fixed load and
domains placed along the z direction. The C++ and Fortran versions are
functionally the same, so the differences between the interfaces can be
checked.

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
