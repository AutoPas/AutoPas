# Gets Eigen, by either (in order priority)
#   1. using an Eigen target a parent project already provides (Eigen3::Eigen or plain Eigen3),
#   2. using an installed/overridden Eigen via find_package
#   3. falls back to the bundled 5.0.1 copy in libs/eigen
# Set Eigen3_ForceBundled=ON to skip 1+2 and always use the bundled copy. If a parent project already provides Eigen 
# this will error, since two copies in one binary would violate the One Definition Rule.
#
# Both `Eigen3::Eigen` (canonical) and `Eigen3` are made to resolve to the same target.

option(Eigen3_ForceBundled "Ignore any provided/installed Eigen and always use the bundled copy." OFF)
mark_as_advanced(Eigen3_ForceBundled)

set(AUTOPAS_EIGEN_MIN_VERSION 5.0.1)

if (NOT Eigen3_ForceBundled AND NOT AUTOPAS_FORCE_ALL_BUNDLED)
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
else ()
    autopas_error_if_forced_bundled_collides(Eigen3 Eigen3_ForceBundled Eigen3::Eigen Eigen3)
endif ()

# Path 3 + fallback: bundled 5.0.1
message(STATUS "Eigen3 - using bundled version 5.0.1")

# Do not add_subdirectory, else we would run configure, build and install. Just define a library from the sources
add_library(
    Eigen3
    OBJECT # this is a header only lib therefore object type is needed
    IMPORTED GLOBAL
)

target_include_directories(Eigen3 SYSTEM INTERFACE "${AUTOPAS_SOURCE_DIR}/libs/eigen")

# Publish the canonical Eigen3::Eigen name alongside the plain Eigen3 target.
autopas_alias_dependency(Eigen3 Eigen3::Eigen)