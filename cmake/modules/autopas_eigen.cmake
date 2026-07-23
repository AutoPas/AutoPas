# Gets Eigen, by either (in order priority)
#   1. using an Eigen target a parent project already provides (Eigen3::Eigen or plain Eigen3),
#   2. using an installed/overridden Eigen via find_package
#   3. falls back to the bundled 3.4.0 copy via FetchContent
# Set Eigen3_ForceBundled=ON to skip 1+2 and always build the bundled copy. This is safe when using AutoPas standalone,
# but dangerous if used within a project with its own Eigen.
#
# Both `Eigen3::Eigen` (canonical) and `Eigen3` are made to resolve to the same target.

option(Eigen3_ForceBundled "Ignore any provided/installed Eigen and always build the bundled copy. Dangerous if used within another project with its own Eigen, safe otherwise." OFF)
mark_as_advanced(Eigen3_ForceBundled)

set(AUTOPAS_EIGEN_MIN_VERSION 3.4.0)

if (NOT Eigen3_ForceBundled)
    # Path 1: reuse an Eigen target a parent project already defined, avoiding a second copy of Eigen in
    # one binary (two copies with differing versions/config macros violate the One Definition Rule (ODR)).
    # We assume the provided target's version, is compatible, as it is hard to check here, and we assume the developer
    # of the parent project will check compatibility. If AutoPas is not compatible, please create an issue.
    if (TARGET Eigen3::Eigen OR TARGET Eigen3)
        message(STATUS "AutoPas: Reusing Eigen provided by parent project")
        autopas_alias_dependency(Eigen3 Eigen3::Eigen)
        return()
    endif ()

    # Path 2: installed or overridden Eigen.
    find_package(Eigen3 ${AUTOPAS_EIGEN_MIN_VERSION} QUIET)
    if (Eigen3_FOUND)
        message(STATUS "Eigen3 - using installed version ${Eigen3_VERSION}")
        autopas_promote_global(Eigen3)
        autopas_promote_global(Eigen3::Eigen)
        autopas_alias_dependency(Eigen3 Eigen3::Eigen)
        return()
    endif ()
    message(STATUS "Eigen3 - no installed version >= ${AUTOPAS_EIGEN_MIN_VERSION} found; using bundled copy")
elseif (TARGET Eigen3::Eigen OR TARGET Eigen3)
    message(
        WARNING
            "Eigen3_ForceBundled=ON but a parent project already provides Eigen; building the bundled copy as well puts two Eigen copies in one binary, risking One Definition Rule (ODR) violations."
    )
endif ()

# Path 3 + fallback: bundled 3.4.0
message(STATUS "Eigen3 - using bundled version 3.4.0 (Release)")

# Enable FetchContent CMake module
include(FetchContent)

# Modern versions of CMake throw warnings from using FetchContent_Populate. We need to still use this to avoid adding
# Eigen's CMakeLists, so we set CMP0169 to OLD to avoid the depreciation warnings.
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0169 OLD)
endif ()

# Build Eigen3 and make the cmake targets available
FetchContent_Declare(
    Eigen3
    URL
        # eigen-master:
        # https://bitbucket.org/eigen/eigen/get/default.zip
        # eigen-3.3.90:
        ${AUTOPAS_SOURCE_DIR}/libs/eigen-3.4.0.zip
    URL_HASH MD5=994092410ba29875184f7725e0371596
)
# In case FetchContent_Populate gets removed:
# Another "hacky" solution is to add the line
#   SOURCE_SUBDIR  _nonexistent_subdir
# to the above FetchContent_Declare, resulting in CMake not finding Eigen's CMakeLists. Then, you can replace all of the
# below (I think) with
#   FetchContent_MakeAvailable(Eigen3)
# This requires at least CMake version 3.18. As of 05/2025, no solution appears to be encouraged, so I hope if
# FetchContent_Populate ever gets removed, there will be some actually encouraged solution but if not, this should work.

# Check if population has already been performed
FetchContent_GetProperties(Eigen3)
if (NOT eigen3_POPULATED) # must be lowercase "eigen3" Fetch the content using previously declared
                          # details
    FetchContent_Populate(Eigen3)

    # Do not add_subdirectory, else we would run configure, build and install Just define a library
    # from the sources
    add_library(
        Eigen3
        OBJECT # this is a header only lib therefore object type is needed
        IMPORTED GLOBAL
    )

    target_include_directories(Eigen3 SYSTEM INTERFACE "${eigen3_SOURCE_DIR}")

    # add_subdirectory(${eigen3_SOURCE_DIR} ${eigen3_BINARY_DIR})
endif ()

# Reset the policy
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0169 NEW)
endif ()

# Publish the canonical Eigen3::Eigen name alongside the plain Eigen3 target.
autopas_alias_dependency(Eigen3 Eigen3::Eigen)