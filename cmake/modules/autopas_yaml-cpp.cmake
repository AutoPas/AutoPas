# Gets yaml-cpp by (in order of priority): reusing a target a parent project provides, an installed
# version via find_package, or the bundled copy in libs/yaml-cpp. Set yaml-cpp_ForceBundled=ON to force bundled.
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
    set(expectedVersion 0.9.0)
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
message(STATUS "yaml-cpp - using bundled version 0.9.0")

# Disable everything we don't need
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Disable yaml-cpp tests")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Disable yaml-cpp contrib")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Disable yaml-cpp tools")
set(YAML_CPP_INSTALL OFF CACHE BOOL "Disable yaml-cpp install rules")

add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/yaml-cpp ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp EXCLUDE_FROM_ALL)

# Hide all remaining internal yaml-cpp cache variables from ccmake/cmake-gui
get_cmake_property(_all_cache_vars CACHE_VARIABLES)
foreach(_var IN LISTS _all_cache_vars)
    if(_var MATCHES "^YAML_")
        mark_as_advanced(${_var})
    endif()
endforeach()

# Disable warnings
target_compile_options(yaml-cpp PRIVATE -w)