.. _usage:

Usage
=====

The library is presented as a single load balancing object ``ALL::ALL``
from ``ALL.hpp``. All classes are encapsulated in the ``ALL`` name space.
Simple (and not so simple) examples are available in the ``/examples`` sub
directory. The simple examples are also documented at length. In cases
where the name `ALL` is misleading, the library is referred to as
`libALL`.

Errors are treated as exceptions that are thrown of a (sub) class of
``ALL::CustomException``. Likely candidates that may throw exceptions are
``balance``,``setup``, ``printVTKoutlines``.

ALL object
----------
The ALL object can be constructed with

.. code-block:: c++

    ALL::ALL<T,W> (
                   const ALL::LB_t method,
                   const int dimension,
                   const T gamma)

Where the type of the boundaries is given by ``T`` and the type of the
work as ``W``, which are usually float or double. The load balancing
method must be one of ``ALL::LB_t::TENSOR``, ``...STAGGERED``,
``...FORCEBASED``, ``...VORONOI``, ``...HISTOGRAM``. Where the Voronoi
method must be enabled at compile time.  There is also a second form where
the initial domain vertices are already provided:

.. code-block:: c++

    ALL::ALL<T,W> (
                   const ALL::LB_t method,
                   const int dimension,
                   const std::vector<ALL::Point<T>> &initialDomain,
                   const T gamma)

Point object
------------
The domain boundaries are given and retrieved as a vector of
``ALL::Point`` s.  These can be constructed via

.. code-block:: c++

    ALL::Point<T> ( const int dimension )
    ALL::Point<T> ( const int dimension, const T *values )
    ALL::Point<T> ( const std::vector<T> &values )

and its elements can be accessed through the ``[]`` operator like a
``std::vector``.

Set up
------
A few additional parameters must be set before the domains can be
balanced.  Which exactly need to be set, is dependent on the method used,
but in general the following are necessary

.. code-block:: c++

    void ALL::ALL<T,W>::setVertices ( std::vector<ALL::Point<T>> &vertices )
    void ALL::ALL<T,W>::setCommunicator ( MPI_Comm comm )
    void ALL::ALL<T,W>::setWork ( const W work )

If the given communicator is not a cartesian communicator the process grid
parameters must also be set beforehand (!) using

.. code-block:: c++

    void ALL::ALL<T,W>::setProcGridParams (
            const std::vector<int> &location,
            const std::vector<int> &size)

with the location of the current process and number of processes in each
direction.

To trace the current domain better an integer tag can be provided, which
is used in the domain output, with

.. code-block:: c++

    void ALL::ALL<T,W>::setProcTag ( int tag )

and an observed minimal domain size in each direction can be set using

.. code-block:: c++

    void ALL::ALL<T,W>::setMinDomainSize (
            const std::vector<T> &minSize )
 
Once these parameters are set call

.. code-block:: c++

    void ALL::ALL<T,W>::setup()

so internal set up routines are run.

Balancing
---------
To create the new domain vertices make sure the current vertices are set
using ``setVertices`` and the work is set using ``setWork`` then call

.. code-block:: c++

    void ALL::ALL<T,W>::balance()

and then the new vertices can be retrieved by using

.. code-block:: c++

    std::vector<ALL::Point<T>> &ALL::ALL<T,W>::getVertices ()

or new neighbors with

.. code-block:: c++

    std::vector<int> &ALL::ALL<T,W>::getNeighbors ()

Special considerations for the Fortran interface
------------------------------------------------
The Fortran interface exists in two versions. Either in the form of

.. code-block:: fortran

    ALL_balance(all_object)

where the all object of type ``type(ALL_t)`` is given as first argument to
each function and as more object oriented functions to the object itself.
So the alternative call would be

.. code-block:: fortran

    all_object%balance()

The function names themselves are similar, i.e. instead of camel case they
are all lowercase except the first ``ALL_`` with underscores between
words. So ``setProcGridParams`` becomes ``ALL_set_proc_grid_parms``. For
more usage info please check the Fortran example programs.

One important difference is the handling of MPI types, which change
between the MPI Fortran 2008 interface as given by the ``mpi_f08`` module
and previous interfaces. In previous interfaces all MPI types are
``INTEGER`` s, but using the ``mpi_f08`` module the communicator is of
``type(MPI_Comm)``. To allow using ``ALL_set_communicator`` directly with
the communicator used in the user's application, build the library with
enabled Fortran 2008 features and this communicator type is expected.

Retrieving information from the balancer is also different, since most
getter return (a reference to) an object itself. The Fortran subroutines
set the values of its arguments. As an example

.. code-block:: c++

    int work = all.getWork();

becomes

.. code-block:: fortran

    integer(c_int) :: work
    call ALL_get_work(work) !or
    !call all%get_work(work)

Since there is no exception handling in Fortran, the error handling from
the libc is borrowed. So there are additional function that retrieve and
reset an error number as well as an additional error message describing
the error. These must be queried and handled by the user. An example of
this usage is shown in the ``ALL_Staggered_f`` example when printing the
VTK outlines.

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
