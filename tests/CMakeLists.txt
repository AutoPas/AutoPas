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

option(AUTOPAS_ENABLE_COVERAGE "Enable code coverage, requires GCC and CMAKE_BUILD_TYPE=Debug." OFF)
message(STATUS "Enable code coverage: " ${AUTOPAS_ENABLE_COVERAGE})

if (AUTOPAS_ENABLE_COVERAGE)
  if (NOT AUTOPAS_BUILD_TESTS OR NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Disabling code coverage.\n CMAKE_BUILD_TYPE=Debug and AUTOPAS_BUILD_TESTS=ON required.")
    set(AUTOPAS_ENABLE_COVERAGE OFF CACHE BOOL "Enable code coverage, requires GCC and CMAKE_BUILD_TYPE=Debug." FORCE)
    message(STATUS "Enable code coverage: " ${AUTOPAS_ENABLE_COVERAGE})
  endif()
  if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "AUTOPAS_ENABLE_COVERAGE only works with GCC. Change your compiler or set AUTOPAS_ENABLE_COVERAGE=OFF")
  endif()
endif()

# includes the gtest modules
include(autopas_gtest)

# add AutoPas tests:
add_subdirectory(testAutopas)
add_subdirectory(MPIParallelAutoPasTests)
