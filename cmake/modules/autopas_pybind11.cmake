set(AUTOPAS_ENABLE_PYTHON_BASED_TUNING
        # Default is OFF just for faster default compilation time.
        OFF
        CACHE
        BOOL "Enables tuning strategies that interface with python. If enabled, pybind11 will be compiled and the json header will be included."
)

if (AUTOPAS_ENABLE_PYTHON_BASED_TUNING)
    # When updating this dependency, update the following expectedVersion and bundledCommit

    set(expectedVersion 2.13.6) # Version of pybind11 which AutoPas is compatible with. This is also the bundled version
    # included with AutoPas. If a non-bundled version is used, it must be compatible with this.
    set(bundledCommit a2e59f0)

    option(pybind11_ForceBundled "Forcibly use the bundled version of pybind11 (v${expectedVersion})" ON)

    if (NOT ${pybind11_ForceBundled})
        find_package(pybind11 ${expectedVersion} QUIET)
        if (pybind11_FOUND)
            message(STATUS "pybind11 - using installed system version ${pybind11_VERSION}")
            # promote target to global visibility to be used elsewhere
            set_target_properties(pybind11::pybind PROPERTIES "IMPORTED_GLOBAL" "TRUE")
            return()
        else ()
            message(STATUS "pybind11 - no system version compatible to version ${expectedVersion} found")
            message(
                    STATUS
                    "pybind11 - if you want to use your version point the cmake variable pybind11_DIR to the directory containing pybind11.cmake in order to pass hints to find_package"
            )
        endif ()
    endif ()

    # system version not found -> install bundled version
    message(STATUS "pybind11 - using bundled version ${expectedVersion} (commit ${bundledCommit})")

    # Enable FetchContent CMake module
    include(FetchContent)

    # Build pybind11 and make the cmake targets available
    FetchContent_Declare(
            pybind11
            URL
            ${AUTOPAS_SOURCE_DIR}/libs/pybind11-${expectedVersion}.zip
            URL_HASH MD5=671deeeaccfccd7c0389e43237c71cf3
            EXCLUDE_FROM_ALL
    )

    FetchContent_MakeAvailable(pybind11)

    # Mark as advanced any options not needed
    MARK_AS_ADVANCED(
            PYBIND11_DISABLE_HANDLE_TYPE_NAME_DEFAULT_IMPLEMENTATION
    )

endif ()


