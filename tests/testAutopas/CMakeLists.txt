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
        $<$<BOOL:${AUTOPAS_ENABLE_COVERAGE}>:gcov>
        autopas
        autopasTools
        molecularDynamicsLibrary
        gmock # gmock includes the gtest target, so we don't need it here.
        Eigen3 # Eigen3 is needed for the tests, but normally not required when using AutoPas.
        $<$<BOOL:${AUTOPAS_ENABLE_HARMONY}>:harmony> # ActiveHarmony is needed for the tests, but normally not required when using AutoPas.
        $<$<BOOL:${AUTOPAS_ENABLE_PYTHON_BASED_TUNING}>:pybind11::embed>
#        $<$<BOOL:${AUTOPAS_ENABLE_TREELITE_BASED_TUNING}>:treelite::treelite>
        treelite::treelite
)


include(GoogleTest)
# queries the compiled executable for tests, this requires the executable to be runnable. if you are
# cross compiling, make sure to properly set CMAKE_CROSSCOMPILING_EMULATOR.
gtest_discover_tests(
        runTests TEST_PREFIX "testAutopas/"
        # increase the discovery timeout for `runTests --gtest_list_tests`
        DISCOVERY_TIMEOUT 60
)
