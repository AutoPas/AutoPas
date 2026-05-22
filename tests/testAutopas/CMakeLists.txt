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
)

# Generate python-based tuning dummy models for tests
if (AUTOPAS_ENABLE_PYTHON_BASED_TUNING)
    set(DT_TEST_MODEL_DIR "${CMAKE_CURRENT_BINARY_DIR}/decisionTreeTuning")
    set(DT_TEST_GENERATE_SCRIPT
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/tuning/tuningStrategy/decisionTreeTuning/generate_test_models.py"
    )

    add_custom_command(TARGET runTests POST_BUILD
        COMMAND "${CMAKE_COMMAND}"
            -DPYTHON="${Python3_EXECUTABLE}"
            -DSCRIPT="${DT_TEST_GENERATE_SCRIPT}"
            -DOUTPUT_DIR="${DT_TEST_MODEL_DIR}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/GenerateDTTestModelsIfMissing.cmake"
        BYPRODUCTS
            "${DT_TEST_MODEL_DIR}/test_model.pkl"
            "${DT_TEST_MODEL_DIR}/test_model_invalid.pkl"
    )

    target_compile_definitions(runTests PRIVATE
        AUTOPAS_DT_TEST_MODEL_DIR="${DT_TEST_MODEL_DIR}"
    )
endif ()

include(GoogleTest)
# queries the compiled executable for tests, this requires the executable to be runnable. if you are
# cross compiling, make sure to properly set CMAKE_CROSSCOMPILING_EMULATOR.
gtest_discover_tests(
        runTests TEST_PREFIX "testAutopas/"
        # Defer --gtest_list_tests to CTest time (not CMake POST_BUILD) to avoid
        # Python interpreter initialisation hanging in the CMake build environment
        # when pybind11::embed is linked.
        DISCOVERY_MODE PRE_TEST
        # increase the discovery timeout for `runTests --gtest_list_tests`
        DISCOVERY_TIMEOUT 60
)
