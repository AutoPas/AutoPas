file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

if (NOT AUTOPAS_INCLUDE_MPI)
    return()
endif()

find_package(MPI)

if (NOT ${MPI_CXX_FOUND})
    message(STATUS "cxx mpi not found, not building runMPITests")
    return()
else ()
    message(STATUS "cxx mpi found: ${MPI_CXX_COMPILER}")
endif ()

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(runMPITests ${MY_SRC})

target_compile_definitions(runMPITests PRIVATE)

target_link_libraries(
    runMPITests
    autopas
    autopasTools
    gmock # gmock includes the gtest target, so we don't need it here.
    Eigen3 # Eigen3 is needed for the tests, but normally not required when using AutoPas.
    harmony # ActiveHarmony is needed for the tests, but normally not required when using AutoPas.
    MPI::MPI_CXX
)

include(GoogleTest)
# queries the compiled executable for tests, this requires the executable to be runnable. if you are
# cross compiling, make sure to properly set CMAKE_CROSSCOMPILING_EMULATOR.
gtest_discover_tests(
    runMPITests TEST_PREFIX "MPIParallelAutoPasTests/"
    # increase the discovery timeout for `runTests --gtest_list_tests`
    DISCOVERY_TIMEOUT 60
)
