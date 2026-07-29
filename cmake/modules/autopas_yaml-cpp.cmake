option(yaml-cpp_ForceBundled "Do not look for an installed version, always use bundled." ON)

if (NOT ${yaml-cpp_ForceBundled})
    set(expectedVersion 0.9.0)
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
