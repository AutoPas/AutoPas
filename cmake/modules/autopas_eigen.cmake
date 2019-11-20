option(Eigen3_ForceBundled "Do not look for an installed version, always use bundled." OFF)

if (NOT ${Eigen3_ForceBundled})
    set(expectedVersion 3.3.4)
    # capital E actually required...
    find_package(Eigen3 ${expectedVersion} QUIET)
    # actually I don't know our minimal supported version but this is the one I tested.
    if (Eigen3_FOUND)
        message(STATUS "Eigen3 - using installed system version ${Eigen3_VERSION}")
        # to later add the alias Eigen3::Eigen needs to be promoted to global visibility
        set_target_properties(Eigen3::Eigen PROPERTIES "IMPORTED_GLOBAL" "TRUE")
        # we need to alias this because aparently make sometimes breaks on '::'
        add_library(Eigen3 ALIAS Eigen3::Eigen)
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
message(STATUS "Eigen3 - using bundled version 3.3.90 (commit 66be6c7)")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract Eigen3
ExternalProject_Add(
    Eigen3_bundled
    URL
        # eigen-master:
        # https://bitbucket.org/eigen/eigen/get/default.zip
        # eigen-3.3.90:
        ${CMAKE_SOURCE_DIR}/libs/eigen-eigen-66be6c76fc01.zip
    URL_HASH MD5=faaf36185ad92b039f7b3f641340dc28
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen-3
    # since we only unpack a header lib src == include
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen-3/include
    # Disable build & install steps
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(Eigen3_bundled source_dir)

add_library(
    Eigen3
    OBJECT      # this is a header only lib therefore object type is needed
    IMPORTED
    GLOBAL
)

add_dependencies(Eigen3 Eigen3_bundled)

# Set libgtest properties
set_target_properties(Eigen3 PROPERTIES "INTERFACE_INCLUDE_DIRECTORIES" "${source_dir}")
