# --------------------------- google benchmark -------------------------------
# Opt-in micro-benchmarks built on top of google/benchmark. Off by default so
# standalone AutoPas builds are not penalized by the extra fetch + compile.

option_conditional(
        AUTOPAS_BUILD_BENCHMARKS
        "Build google-benchmark based micro-benchmarks (will be downloaded if not available)"
        OFF
)

if (NOT AUTOPAS_BUILD_BENCHMARKS)
    return()
endif ()

# Fetch benchmark::benchmark from the bundled zip under libs/.
include(autopas_benchmark)

# Each micro-benchmark lives in a dedicated file to keep compile times low.
add_executable(LJFunctorHWYBench LJFunctorHWYBench.cpp)

target_link_libraries(
        LJFunctorHWYBench
        PRIVATE
        autopas
        molecularDynamicsLibrary
        benchmark::benchmark
)

# Debug-only tool for the sorting-threshold calibration investigation (debug/sorting-threshold-
# cache-residency branch). Does not depend on google-benchmark. Not part of the thesis-cited
# benchmark dataset.
add_executable(SortingCacheResidencyDebug SortingCacheResidencyDebug.cpp)

target_link_libraries(
        SortingCacheResidencyDebug
        PRIVATE
        autopas
        molecularDynamicsLibrary
)

# MoleculeLJ lives in the single-site MD mode. The tests rely on the same
# MD_FLEXIBLE_MODE=0 setting inherited from molecularDynamicsLibrary.
