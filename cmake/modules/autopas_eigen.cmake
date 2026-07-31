option(Eigen3_ForceBundled "Do not look for an installed version, always use bundled." ON)

if (NOT ${Eigen3_ForceBundled})
    set(expectedVersion 3.4.0)
    # capital E actually required...
    find_package(Eigen3 ${expectedVersion} QUIET)
    # actually I don't know our minimal supported version but this is the one I tested.
    if (Eigen3_FOUND)
        message(STATUS "Eigen3 - using installed system version ${Eigen3_VERSION}")
        # to later add the alias Eigen3::Eigen needs to be promoted to global visibility
        set_target_properties(Eigen3::Eigen PROPERTIES "IMPORTED_GLOBAL" "TRUE")
        # we need to alias this because aparently make sometimes breaks on '::'
        add_library(Eigen3 ALIAS Eigen3::Eigen)
        # return, as we have found the target
        return()
    else ()
        message(STATUS "Eigen3 - no system version compatible to version ${expectedVersion} found")
        message(
            STATUS
                "Eigen3 - if you want to use your version point the cmake variable Eigen3_DIR to the directory containing Eigen3Config.cmake in order to pass hints to find_package"
        )
    endif ()
endif ()

# system version not found -> install bundled version
message(STATUS "Eigen3 - using bundled version 3.4.0 (Release)")

# Do not add_subdirectory, else we would run configure, build and install. Just define a library from the sources
add_library(
    Eigen3
    OBJECT # this is a header only lib therefore object type is needed
    IMPORTED GLOBAL
)

target_include_directories(Eigen3 SYSTEM INTERFACE "${AUTOPAS_SOURCE_DIR}/libs/eigen")
