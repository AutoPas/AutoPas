option(simde_ForceBundled "Do not look for an installed version, always use bundled." ON)

if (NOT ${simde_ForceBundled})
    set(expectedVersion 0.8.0)
    # capital E actually required...
    find_package(simde ${expectedVersion} QUIET)
    # actually I don't know our minimal supported version but this is the one I tested.
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
message(STATUS "simde - using bundled version 0.8.0 (commit e651ec3)")

# Enable FetchContent CMake module
include(FetchContent)

# Build Eigen3 and make the cmake targets available
FetchContent_Declare(
    simde
    URL
        ${AUTOPAS_SOURCE_DIR}/libs/simde-0.8.0-rc1.zip
    URL_HASH MD5=d95da9e9c5377e96128ea87caa9be30d
)

# Check if population has already been performed
FetchContent_GetProperties(simde)
if (NOT simde_POPULATED) # must be lowercase "eigen3" Fetch the content using previously declared
                          # details
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
