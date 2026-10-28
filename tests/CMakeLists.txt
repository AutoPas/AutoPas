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

# Shared testing helpers
add_library(testingHelpers INTERFACE)
target_include_directories(testingHelpers INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/testingHelpers"
)

# add AutoPas tests:
add_subdirectory(testAutopas)
add_subdirectory(MPIParallelAutoPasTests)
