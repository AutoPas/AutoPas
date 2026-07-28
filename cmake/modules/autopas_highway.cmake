# Disable unnecessary Highway features that bloat the build and cache
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Disable Highway tests" FORCE) # Force OFF to avoid conflicting test targets
set(HWY_ENABLE_CONTRIB OFF CACHE BOOL "Disable Highway contrib")
set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "Disable Highway examples")
set(HWY_ENABLE_INSTALL OFF CACHE BOOL "Disable Highway install")

# Hide all Highway-specific CMake options from the standard ccmake view
mark_as_advanced(
        HWY_CMAKE_ARM7           # Manually forces ARMv7 build (usually auto-detected)
        HWY_CMAKE_HEADER_ONLY    # Forces Highway to not compile its .cc files
        HWY_CMAKE_RVV            # Manually enables RISC-V Vector extensions
        HWY_CMAKE_SSE2           # Manually enables SSE2 for x86
        HWY_DISABLE_FUTEX        # Disables fast-mutexes in Highway's thread pool
        HWY_ENABLE_CONTRIB       # Builds extra Highway algorithms
        HWY_ENABLE_EXAMPLES      # Builds Highway examples
        HWY_ENABLE_INSTALL       # Installs Highway globally during `make install`
        HWY_ENABLE_TESTS         # Builds Highway's own test suite (forced OFF above)
        HWY_TEST_STANDALONE      # Builds Highway tests without Google Test
        HWY_WARNINGS_ARE_ERRORS  # Highway developer flag to treat compiler warnings as errors
)

add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/highway ${CMAKE_CURRENT_BINARY_DIR}/highway EXCLUDE_FROM_ALL)
