option(yaml-cpp_ForceBundled "Do not look for an installed version, always use bundled." ON)

if (NOT ${yaml-cpp_ForceBundled})
    set(expectedVersion 0.8.0)
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
# This is not a stable release, (is something after 0.8.0), but is required for compatibility with CMake 4.0.
message(STATUS "yaml-cpp - using bundled version 2f86d13")

# Enable FetchContent CMake module
include(FetchContent)

# Build yaml-cpp and make the cmake targets available
FetchContent_Declare(
    yaml-cpp
    URL
        # yaml-cpp-master:
        # https://github.com/jbeder/yaml-cpp/archive/refs/heads/master.zip
        # commit 2f86d13:
        ${AUTOPAS_SOURCE_DIR}/libs/yaml-cpp-2f86d13.zip
    URL_HASH MD5=d402b60e57c14fcb30138c5b28a333d1
    # needed to compile with ninja
)

# Disable everything we don't need
option(YAML_CPP_BUILD_TESTS "" OFF)
option(YAML_CPP_BUILD_CONTRIB "" OFF)
option(YAML_CPP_BUILD_TOOLS "" OFF)

# hide options from ccmake
mark_as_advanced(
    YAML_BUILD_SHARED_LIBS
    YAML_CPP_BUILD_CONTRIB
    YAML_CPP_BUILD_TESTS
    YAML_CPP_BUILD_TOOLS
    YAML_CPP_CLANG_FORMAT_EXE
)

FetchContent_MakeAvailable(yaml-cpp)

# Disable warnings
target_compile_options(yaml-cpp PRIVATE -w)
