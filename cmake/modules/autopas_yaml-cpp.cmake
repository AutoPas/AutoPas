# Gets yaml-cpp by (in order of priority): reusing a target a parent project provides, an installed
# version via find_package, or the bundled copy. Set yaml-cpp_ForceBundled=ON to force bundled.
option(yaml-cpp_ForceBundled "Ignore any provided/installed yaml-cpp and always use the bundled copy." OFF)
mark_as_advanced(yaml-cpp_ForceBundled)

if (NOT yaml-cpp_ForceBundled)
    # Path 1: reuse a yaml-cpp target a parent project already defined.
    if (TARGET yaml-cpp OR TARGET yaml-cpp::yaml-cpp)
        message(STATUS "AutoPas: Reusing yaml-cpp provided by parent project")
        autopas_alias_dependency(yaml-cpp yaml-cpp::yaml-cpp)
        return()
    endif ()
    # Path 2: installed version; the version arg enforces our minimum. Use whatever compatible version
    # find_package returns - it creates a target either way, so rejecting it here would leave a
    # `yaml-cpp::yaml-cpp` that collides with the bundled copy's own alias below.
    set(expectedVersion 0.8.0)
    find_package(yaml-cpp ${expectedVersion} QUIET)
    if (yaml-cpp_FOUND)
        message(STATUS "yaml-cpp - using installed system version ${yaml-cpp_VERSION}")
        autopas_promote_global(yaml-cpp)
        autopas_promote_global(yaml-cpp::yaml-cpp)
        autopas_alias_dependency(yaml-cpp yaml-cpp::yaml-cpp)
        return()
    endif ()
    message(STATUS "yaml-cpp - no system version >= ${expectedVersion} found; using bundled copy")
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
