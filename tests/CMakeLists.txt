# make test optional. Useful also when no network connection is available
option_conditional(
    AUTOPAS_BUILD_TESTS
    "Build tests based on the Google Test Framework (will be downloaded if not available)"
    AUTOPAS_STANDALONE_BUILD
)

# do nothing if nothing should be done
if (NOT AUTOPAS_BUILD_TESTS)
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building tests.")
    return()
endif ()

# includes the gtest modules
include(autopas_gtest)

# Find and include Python headers and libraries
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Include the Python headers
include_directories(${Python3_INCLUDE_DIRS})

# add AutoPas tests:
add_subdirectory(testAutopas)
add_subdirectory(MPIParallelAutoPasTests)
