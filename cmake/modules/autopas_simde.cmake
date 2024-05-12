option(simde_ForceBundled "Do not look for an installed version, always use bundled." ON)

if (NOT ${simde_ForceBundled})
    set(expectedVersion 0.8.2)
    find_package(simde ${expectedVersion} QUIET)
    if (simde_FOUND)
        message(STATUS "simde - using installed system version ${simde_VERSION}")
        # return, as we have found the target
        return()
    else ()
        message(STATUS "simde - no system version compatible to version ${expectedVersion} found")
        message(
            STATUS
                "simde - if you want to use your version point the cmake variable simde_DIR to the directory containing simdeConfig.cmake in order to pass hints to find_package"
        )
    endif ()
endif ()

# system version not found -> install bundled version
message(STATUS "simde - using bundled version 0.8.2 (commit 71fd833)")

# Enable FetchContent CMake module
include(FetchContent)

# Build SIMDe and make the cmake targets available
FetchContent_Declare(
    simde
    URL
        ${AUTOPAS_SOURCE_DIR}/libs/simde-0.8.2.zip
    URL_HASH MD5=d0dadc0fbc56bdbf33d7eaf28cf7e439
)

# Check if population has already been performed
FetchContent_GetProperties(simde)
if (NOT simde_POPULATED)
    FetchContent_Populate(simde)

    # Do not add_subdirectory, else we would run configure, build and install Just define a library
    # from the sources
    add_library(
        simde
        OBJECT # this is a header only lib therefore object type is needed
        IMPORTED GLOBAL
    )

    target_include_directories(simde SYSTEM INTERFACE "${simde_SOURCE_DIR}")

    # add_subdirectory(${eigen3_SOURCE_DIR} ${eigen3_BINARY_DIR})
endif ()
