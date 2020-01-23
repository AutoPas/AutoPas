option(yaml-cpp_ForceBundled "Do not look for an installed version, always use bundled." OFF)

if (NOT ${yaml-cpp_ForceBundled})
    set(expectedVersion 0.5.2)
    # first try: check if we find any installed version
    find_package(yaml-cpp ${expectedVersion} QUIET)
    if (yaml-cpp_FOUND)
        message(STATUS "yaml-cpp - using installed system version ${yaml-cpp_VERSION}")
        # return here, as we have now found and imported the target.
        return()
    else ()
        message(
            STATUS "yaml-cpp - no system version compatible to version ${expectedVersion} found"
        )
        message(
            STATUS
                "yaml-cpp - if you want to use your version point the cmake variable yaml-cpp_DIR to the directory containing  yaml-cpp-config.cmake in order to pass hints to find_package"
        )
    endif ()
endif ()

# system version not found -> install bundled version
message(STATUS "yaml-cpp - using bundled version 0.6.3 (commit 72f699f)")

# Enable FetchContent CMake module
include(FetchContent)

# Build spdlog and make the cmake targets available
FetchContent_Declare(
    yaml-cpp
    URL
    # yaml-cpp-master:
    # https://github.com/jbeder/yaml-cpp/archive/master.zip
    # commit 72f699f:
    ${AUTOPAS_SOURCE_DIR}/libs/yaml-cpp-master.zip
    URL_HASH MD5=6186f7ea92b907e9128bc74c96c1f791
    # needed to compile with ninja
)

# Disable everything we don't need and set build type to release. Also disable warnings.
set(YAML_CPP_BUILD_TESTS CACHE BOOL OFF)
set(YAML_CPP_BUILD_CONTRIB CACHE BOOL OFF)
set(YAML_CPP_BUILD_TOOLS CACHE BOOL OFF)

# Things that were previously set, but we can't set now because the library is no longer
# a sub-build and we can't pass seperate CMake args
# Maybe this will be possible with this feature from CMake 3.17.0
# https://gitlab.kitware.com/cmake/cmake/issues/19854
# BUILD_TESTING=OFF
# CMAKE_BUILD_TYPE=RELEASE
# CMAKE_CXX_FLAGS=-w

FetchContent_MakeAvailable(yaml-cpp)