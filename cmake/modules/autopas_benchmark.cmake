# Gets Google Benchmark by (in order of priority): reusing a `benchmark` target a parent project provides, an installed
# version via find_package, or the bundled copy in libs/benchmark. Set benchmark_ForceBundled=ON to force bundled.
option(benchmark_ForceBundled "Ignore any provided/installed Google Benchmark and always use the bundled copy." OFF)
mark_as_advanced(benchmark_ForceBundled)

if (NOT benchmark_ForceBundled AND NOT AUTOPAS_FORCE_ALL_BUNDLED)
    set(expectedVersion 1.9.5)
    # Path 1: reuse a Google Benchmark target a parent project already defined.
    if (TARGET benchmark OR TARGET benchmark::benchmark)
        message(STATUS "AutoPas: Reusing Google Benchmark provided by parent project")
        autopas_warn_if_parent_version_too_old(Benchmark benchmark ${expectedVersion} benchmark benchmark::benchmark)
        autopas_alias_dependency(benchmark benchmark::benchmark)
        return()
    endif ()
    # Path 2: installed version; the version arg enforces our minimum.
    find_package(benchmark ${expectedVersion} QUIET)
    if (benchmark_FOUND)
        message(STATUS "Google Benchmark - using installed version ${benchmark_VERSION}")
        autopas_promote_global(benchmark)
        autopas_promote_global(benchmark::benchmark)
        autopas_alias_dependency(benchmark benchmark::benchmark)
        return()
    endif ()
    message(STATUS "Google Benchmark - no installed version >= ${expectedVersion} found; using bundled copy")
else ()
    autopas_error_if_forced_bundled_collides(Benchmark benchmark_ForceBundled benchmark benchmark::benchmark)
endif ()

# Path 3 + fallback: bundled version.
message(STATUS "Google Benchmark - using bundled version 1.9.5")

find_package(Threads REQUIRED)

# Disable unnecessary Google Benchmark features that bloat the build and cache.
# Tests must stay off, otherwise benchmark creates a gtest target in conflict with our own.
set(BENCHMARK_ENABLE_TESTING        OFF CACHE INTERNAL "Disable benchmark tests")
set(BENCHMARK_ENABLE_GTEST_TESTS    OFF CACHE INTERNAL "Disable benchmark gtest unit tests")
set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF CACHE INTERNAL "Disable benchmark assembly verification tests")
set(BENCHMARK_ENABLE_DOXYGEN        OFF CACHE INTERNAL "Disable benchmark Doxygen documentation")
set(BENCHMARK_ENABLE_WERROR         OFF CACHE INTERNAL "Disable building release with -Werror")
set(BENCHMARK_FORCE_WERROR          OFF CACHE INTERNAL "Disable forcing -Werror regardless of compiler")
set(BENCHMARK_ENABLE_INSTALL        OFF CACHE INTERNAL "Disable benchmark installation rules")
set(BENCHMARK_INSTALL_DOCS          OFF CACHE INTERNAL "Disable benchmark documentation installation")
set(BENCHMARK_INSTALL_TOOLS         OFF CACHE INTERNAL "Disable benchmark Python tools installation")
set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF CACHE INTERNAL "Disable downloading unmet dependencies time")

mark_as_advanced(
        BENCHMARK_BUILD_32_BITS      # Builds 32-bit version of the benchmark library
        BENCHMARK_ENABLE_EXCEPTIONS  # Enables use of C++ exceptions in benchmark library
        BENCHMARK_ENABLE_LIBPFM      # Enables performance counters provided by libpfm
        BENCHMARK_ENABLE_LTO         # Enables Link Time Optimization (LTO)
        BENCHMARK_USE_BUNDLED_GTEST  # Uses bundled GoogleTest instead of find_package
        BENCHMARK_USE_LIBCXX         # Builds and tests using LLVM libc++ standard library
)

# This is a small workaround to prevent Google Benchmark from printing the AutPas version instead of its own
set(__get_git_version INCLUDED)
function(get_git_version var)
    set(${var} "v0.0.0" PARENT_SCOPE)
endfunction()

add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/benchmark ${CMAKE_CURRENT_BINARY_DIR}/benchmark
        EXCLUDE_FROM_ALL)

# expose the namespaced name too, so consumers can link either benchmark or benchmark::benchmark
autopas_alias_dependency(benchmark benchmark::benchmark)