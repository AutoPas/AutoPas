# The VTK error test case, where the output directory is not accessible.
include(FindUnixCommands)
set(OUT_DIR vtk_outline)
add_test(
    NAME fortran_error_vtk_setup
    COMMAND ${BASH} -c "mkdir -p ${OUT_DIR} && chmod -x ${OUT_DIR}; exit $?"
    )
set_tests_properties(fortran_error_vtk_setup PROPERTIES
    LABELS "ALL;fortran_error"
    FIXTURES_SETUP fortran_error_filesystem
    )
add_test(
    NAME fortran_error_vtk_cleanup
    COMMAND ${BASH} -c "chmod +x ${OUT_DIR} && rm -rf ${OUT_DIR}"
    )
set_tests_properties(fortran_error_vtk_cleanup PROPERTIES
    LABELS "ALL;fortran_error"
    FIXTURES_CLEANUP fortran_error_filesystem
    )
add_test(
    NAME fortran_error_vtk
    COMMAND ${BASH} -c "mpirun -n 1 $<TARGET_FILE:ALL_Staggered_f>; exit $?"
    )
set_tests_properties(fortran_error_vtk PROPERTIES
    LABELS "ALL;fortran_error"
    FIXTURES_REQUIRED fortran_error_filesystem
    WILL_FAIL TRUE
    )


# vim: et sw=4 ts=4
