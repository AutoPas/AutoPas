# Gets Highway by (in order of priority): reusing a `hwy` target a parent project provides, an installed
# version via find_package, or the bundled copy in libs/highway. Set highway_ForceBundled=ON to force bundled.
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
    set(expectedVersion 1.4.0)
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
message(STATUS "Highway - using bundled version 1.4.0")

# Disable unnecessary Highway features that bloat the build and cache.
# Tests must stay off, otherwise Highway creates a gtest target in conflict with our own.
set(HWY_ENABLE_TESTS OFF CACHE INTERNAL "Disable Highway tests")
set(HWY_ENABLE_CONTRIB OFF CACHE INTERNAL "Disable Highway contrib")
set(HWY_ENABLE_EXAMPLES OFF CACHE INTERNAL "Disable Highway examples")
set(HWY_ENABLE_INSTALL OFF CACHE INTERNAL "Disable Highway install")

# Hide all Highway-specific CMake options from the standard ccmake view
mark_as_advanced(
        HWY_CMAKE_ARM7           # Manually forces ARMv7 build (usually auto-detected)
        HWY_CMAKE_HEADER_ONLY    # Forces Highway to not compile its .cc files
        HWY_CMAKE_RVV            # Manually enables RISC-V Vector extensions
        HWY_CMAKE_SSE2           # Manually enables SSE2 for x86
        HWY_DISABLE_FUTEX        # Disables fast-mutexes in Highway's thread pool
        HWY_TEST_STANDALONE      # Builds Highway tests without Google Test
        HWY_WARNINGS_ARE_ERRORS  # Highway developer flag to treat compiler warnings as errors
)

add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/highway ${CMAKE_CURRENT_BINARY_DIR}/highway EXCLUDE_FROM_ALL)

# expose the namespaced name too, so consumers can link either hwy or hwy::hwy
autopas_alias_dependency(hwy hwy::hwy)