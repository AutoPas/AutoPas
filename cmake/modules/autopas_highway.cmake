# Gets Highway by (in order of priority): reusing a `hwy` target a parent project provides, an installed
# version via find_package, or the bundled 1.3.0 copy. Set highway_ForceBundled=ON to force bundled.
option(highway_ForceBundled "Ignore any provided/installed Highway and always use the bundled copy." OFF)
mark_as_advanced(highway_ForceBundled)

if (NOT highway_ForceBundled)
    # Path 1: reuse a Highway target a parent project already defined.
    if (TARGET hwy OR TARGET hwy::hwy)
        message(STATUS "AutoPas: Reusing Google Highway provided by parent project")
        autopas_alias_dependency(hwy hwy::hwy)
        return()
    endif ()
    # Path 2: installed version; the version arg enforces our minimum.
    set(expectedVersion 1.3.0)
    find_package(hwy ${expectedVersion} QUIET)
    if (hwy_FOUND)
        message(STATUS "Highway - using installed version ${hwy_VERSION}")
        autopas_promote_global(hwy)
        autopas_promote_global(hwy::hwy)
        autopas_alias_dependency(hwy hwy::hwy)
        return()
    endif ()
    message(STATUS "Highway - no installed version >= ${expectedVersion} found; using bundled copy")
endif ()

# Path 3 + fallback: bundled version.
message(STATUS "Highway - using bundled version 1.3.0")
include(FetchContent)

FetchContent_Declare(
    autopas_highway
    URL ${AUTOPAS_SOURCE_DIR}/libs/highway-1.3.0.zip
)

# We must force this to ensure that Highway does not create a gtest target in conflict with our own.
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Disable Highway tests" FORCE)
FetchContent_MakeAvailable(autopas_highway)

# expose the namespaced name too, so consumers can link either hwy or hwy::hwy
autopas_alias_dependency(hwy hwy::hwy)