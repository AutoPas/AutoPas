if (NOT AUTOPAS_INCLUDE_MPI)
    message(STATUS "runMPITests - MPI not included. Skipping!")
    return()
endif()
# at this point MPI should already have been found via the autopas_mpi module

message(STATUS "runMPITests - adding tests.")

# all sources in this folder and sub folders are considered relevant for this target
file(
        GLOB_RECURSE
        MY_SRC
        "*.cpp"
        "*.h"
)
add_executable(runMPITests ${MY_SRC})

target_link_libraries(
    runMPITests
    autopas
    autopasTools
    gmock
    Eigen3 # Eigen3 is needed for the tests, but normally not required when using AutoPas.
    harmony # ActiveHarmony is needed for the tests, but normally not required when using AutoPas.
    MPI::MPI_CXX
)
set(AUTOPAS_TESTS_NUM_MPI_RANKS 3 CACHE STRING "Number of MPI Ranks to run AutoPas tests with.")
set_property(TARGET runMPITests PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${AUTOPAS_TESTS_NUM_MPI_RANKS}')

include(GoogleTest)
# queries the compiled executable for tests, this requires the executable to be runnable. if you are
# cross compiling, make sure to properly set CMAKE_CROSSCOMPILING_EMULATOR.
gtest_discover_tests(
    runMPITests TEST_PREFIX "MPIParallelAutoPasTests/"
    # increase the discovery timeout for `runTests --gtest_list_tests`
    DISCOVERY_TIMEOUT 60
)
