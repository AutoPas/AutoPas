file(
        GLOB_RECURSE
        MY_SRC
        "*.cpp"
        "*.h"
)

add_executable(runTests ${MY_SRC})

target_compile_definitions(runTests PRIVATE)

target_link_libraries(
        runTests
        PUBLIC
        $<$<BOOL:${AUTOPAS_ENABLE_COVERAGE}>:gcov>
        autopas
        generators
        molecularDynamicsLibrary
        gmock # gmock includes the gtest target, so we don't need it here.
        Eigen3 # Eigen3 is needed for the tests, but normally not required when using AutoPas.
        $<$<BOOL:${AUTOPAS_ENABLE_HARMONY}>:harmony> # ActiveHarmony is needed for the tests, but normally not required when using AutoPas.
)

# Kokkos (and its CUDA headers) do not build against libc++, so Clang-based Kokkos builds have to
# link explicitly against libstdc++ and libm.
if (AUTOPAS_ENABLE_KOKKOS AND (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
    target_link_libraries(runTests PRIVATE stdc++ m)
endif ()

include(GoogleTest)
# queries the compiled executable for tests, this requires the executable to be runnable. if you are
# cross compiling, make sure to properly set CMAKE_CROSSCOMPILING_EMULATOR.
gtest_discover_tests(
        runTests TEST_PREFIX "testAutopas/"
        # Use PRE_TEST discovery mode to avoid running MPI executable during build
        # This prevents SIGPIPE errors in CI when MPI and sanitizers are enabled (this is needed still, even with these
        # non-MPI tests)
        DISCOVERY_MODE PRE_TEST
        # increase the discovery timeout for `runTests --gtest_list_tests`
        DISCOVERY_TIMEOUT 60
)
