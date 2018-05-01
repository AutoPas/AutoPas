# make test optional. Useful also when no network connection is available
option(BUILD_TESTS "Build tests based on the Google Test Framework (will be downloaded if not available)" ON)

# do nothing if nothing should be done
if (NOT BUILD_TESTS)
    return()
endif ()

#includes the gtest modules
include(autopas_gtest)

# add AutoPas tests:
add_subdirectory(testAutopas)