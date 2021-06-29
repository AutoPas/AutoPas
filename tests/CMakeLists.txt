# make test optional. Useful also when no network connection is available
cmake_dependent_option(
    AUTOPAS_BUILD_TESTS
    "Build tests based on the Google Test Framework (will be downloaded if not available)"
    ON
    AUTOPAS_STANDALONE_BUILD
    OFF
)

# do nothing if nothing should be done
if (NOT AUTOPAS_BUILD_TESTS)
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building tests.")
    return()
endif ()

# includes the gtest modules
include(autopas_gtest)

# add AutoPas tests:
add_subdirectory(testAutopas)
add_subdirectory(MPIParallelAutoPasTests)
