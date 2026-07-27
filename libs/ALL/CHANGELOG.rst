=========
Changelog
=========

Version 1.0
-----------

Version 1.0.0
*************

Version 0.9
-----------

Version 0.9.3
*************

- Bug: ALL::Point.dist_plane(..) had a type bug. This manifestes with a cuda
  compiler even without explicit usage of the function. Fixed and added
  corresponding test.
- Add support for external voro++
- Bug: fix race condition in parallel test execution
- Allow custom install locations in CMake via GNUINstallDirs and add so version
- Bug: removed unncessary recreation of internal communicators in STAGGERED and tensor
- Bug: offset minimum sizes are now considered in HISTOGRAM
- Feature: function to compute load efficiency
- Bug: Fixed wrong behavior when computing floor values of negative values
- several other small fixes and corrections

Version 0.9.2
*************

- *CHANGE*: The old CMake policy ``CMP0004``, which was required for some MPI
  installations, is no longer set. If this is still required: Please contact
  us.
- *CHANGE*: Fortran modules are now installed in ``/lib``, instead of
  ``/include/modules``.
- Feature: Integration tests are only generated if ``CM_ALL_TESTS_INTEGRATION``
  is set.
- Feature: Example CMake and Make projects for integrating ALL into the build
  process.
- Bug: CMake dependencies between targets and link and include inheritance
  corrected.
- Feature: Add ALL_Tensor feature test.
- Bug: VTK fails to include <limits>, so we take care of that (for GCC 11 and
  VTK 9.3.0)
- Bug: setVertices was required for printVTKoutlines to output updated vertices.
- Feature: Fortran API supports retrieval of error string
- Feature: Testing against VTK 7 as well as 9
- Bug: A nested namespace definition was removed, which is only allowed in C++17.

Version 0.9.1
*************

- Bug: Wrong gamma calculation for staggered grid.
- Bug: No optimal gamma calculation for tensor method.

Version 0.9.0
*************
Initial release
