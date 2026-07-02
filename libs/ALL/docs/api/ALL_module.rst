.. _ALLModule:

ALL Fortran module
==================

Due to the case insensitive nature of Fortran and the used toolchain
(doxygen and sphinx), all variables and function names are lowercase on
this page. We do, however, use only upper case for global constants, such
as the load balancing methods (``ALL_STAGGERED`` etc.) and the function
names are usually written similarly to MPI as ``ALL_get_work``. The object
is of ``type(ALL_t)``.

Methods
-------
Available methods are:

.. doxygenvariable:: all_forcebased
.. doxygenvariable:: all_histogram
.. doxygenvariable:: all_staggered
.. doxygenvariable:: all_tensor
.. doxygenvariable:: all_voronoi


Errors
------
The exceptions map to the following error constants:

.. doxygenvariable:: all_error_filesystem
.. doxygenvariable:: all_error_generic
.. doxygenvariable:: all_error_internal
.. doxygenvariable:: all_error_invalidargument
.. doxygenvariable:: all_error_invalidcommtype
.. doxygenvariable:: all_error_outofbounds
.. doxygenvariable:: all_error_pointdimensionmissmatch

And the error message is returned in a string of length
``ALL_ERROR_LENGTH``.

Functions
---------
These function all have a corresponding call in the object, with the
``ALL_`` stripped from the function name. So with an object named
``balancer`` the function ``ALL_get_work(balancer, work)`` is identical to
``balancer%get_work(work)``. The following functions are available:

.. doxygenfunction:: all_balance
.. doxygenfunction:: all_finalize
.. doxygenfunction:: all_setup

Getters
*******
.. doxygenfunction:: all_get_dimension
.. doxygenfunction:: all_get_gamma
.. doxygenfunction:: all_get_length_of_work
.. doxygenfunction:: all_get_neighbors
.. doxygenfunction:: all_get_number_of_neighbors
.. doxygenfunction:: all_get_number_of_vertices
.. doxygenfunction:: all_get_prev_vertices
.. doxygenfunction:: all_get_vertices
.. doxygenfunction:: all_get_vertices_alloc
.. doxygenfunction:: all_get_work
.. doxygenfunction:: all_get_work_array

Setters
*******
.. doxygenfunction:: all_set_gamma
.. doxygenfunction:: all_set_method_data_histgram
.. doxygenfunction:: all_set_min_domain_size
.. doxygenfunction:: all_set_proc_grid_params
.. doxygenfunction:: all_set_proc_tag
.. doxygenfunction:: all_set_sys_size
.. doxygenfunction:: all_set_vertices
.. doxygenfunction:: all_set_work
.. doxygenfunction:: all_set_work_multi

Output
******
.. doxygenfunction:: all_print_vtk_outlines
.. doxygenfunction:: all_print_vtk_vertices

Error handling
**************
.. doxygenfunction:: all_error
.. doxygenfunction:: all_error_description
.. doxygenfunction:: all_reset_error

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
