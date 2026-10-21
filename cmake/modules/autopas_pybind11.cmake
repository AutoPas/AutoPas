set(AUTOPAS_ENABLE_PYTHON_BASED_TUNING
        # Default is OFF just for faster default compilation time.
        OFF
        CACHE
        BOOL "Enables tuning strategies that interface with python. If enabled, pybind11 will be compiled and the json header will be included."
)

if (AUTOPAS_ENABLE_PYTHON_BASED_TUNING)
    # Gets pybind11 by (in order of priority): reusing a target a parent project provides, an installed
    # version via find_package, or the bundled copy in libs/pybind11. Set pybind11_ForceBundled=ON to force bundled.

    # When updating this dependency, update the following expectedVersion and bundledCommit

    set(expectedVersion 3.1.0) # Version of pybind11 which AutoPas is compatible with. This is also the bundled version
    # included with AutoPas. If a non-bundled version is used, it must be compatible with this.
    set(bundledCommit 97bf890)

    option(pybind11_ForceBundled "Ignore any provided/installed pybind11 and always use the bundled copy." OFF)
    mark_as_advanced(pybind11_ForceBundled)

    set(PYBIND11_FINDPYTHON NEW)

    if (NOT pybind11_ForceBundled AND NOT AUTOPAS_FORCE_ALL_BUNDLED)
        # Path 1: reuse a pybind11 a parent project already defined.
        if (TARGET pybind11::embed)
            message(STATUS "AutoPas: Reusing pybind11 provided by parent project")
            autopas_warn_if_parent_version_too_old(pybind11 pybind11 ${expectedVersion})
            return()
        endif ()
        # Path 2: installed version
        find_package(pybind11 ${expectedVersion} QUIET)
        if (pybind11_FOUND)
            message(STATUS "pybind11 - using installed system version ${pybind11_VERSION}")
            # promote the targets AutoPas links to global visibility to be used elsewhere
            autopas_promote_global(pybind11::embed)
            return()
        endif ()
        message(STATUS "pybind11 - no system version >= ${expectedVersion} found; using bundled copy")
    else ()
        autopas_error_if_forced_bundled_collides(pybind11 pybind11_ForceBundled pybind11::embed)
    endif ()

    # Path 3 + fallback: bundled version.
    message(STATUS "pybind11 - using bundled version ${expectedVersion} (commit ${bundledCommit})")

    # Build pybind11 and make the cmake targets available
    add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/pybind11 EXCLUDE_FROM_ALL)

    # Hide all remaining internal pybind11 cache variables from ccmake/cmake-gui
    get_cmake_property(_all_cache_vars CACHE_VARIABLES)
    foreach(_var IN LISTS _all_cache_vars)
        if(_var MATCHES "^PYBIND11_")
            mark_as_advanced(${_var})
        endif()
    endforeach()

endif ()