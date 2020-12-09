file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)
if (AUTOPAS_ENABLE_CUDA)
    file(
        GLOB_RECURSE
        CU_SRC
        "*.cu"
        "*.cuh"
    )
endif ()

option(AUTOPAS_USE_TEST_GLOB "Testing: Use one cpp file that includes all other cpp files. Reduces the amount of needed template instantiations." OFF)

if(AUTOPAS_USE_TEST_GLOB)
    set(NEW_LIST ${MY_SRC})
    list(TRANSFORM NEW_LIST PREPEND "#include \"")
    list(TRANSFORM NEW_LIST APPEND "\"\n")
    set(AUTOPAS_TEST_FILES ${CMAKE_CURRENT_BINARY_DIR}/test_glob.cpp)
    file(WRITE ${AUTOPAS_TEST_FILES} ${NEW_LIST})
else()
    set(AUTOPAS_TEST_FILES ${MY_SRC})
endif()

add_executable(runTests $<$<BOOL:${AUTOPAS_ENABLE_CUDA}>:${CU_SRC}> ${AUTOPAS_TEST_FILES})

target_compile_definitions(runTests PRIVATE)

target_link_libraries(
    runTests
    autopas
    autopasTools
    gmock # gmock includes the gtest target, so we don't need it here.
    Eigen3 # Eigen3 is needed for the tests, but normally not required when using AutoPas.
    harmony # ActiveHarmony is needed for the tests, but normally not required when using AutoPas.
)

include(GoogleTest)
# queries the compiled executable for tests, this requires the executable to be runnable. if you are
# cross compiling, make sure to properly set CMAKE_CROSSCOMPILING_EMULATOR.
gtest_discover_tests(
    runTests TEST_PREFIX "testAutopas/"
    # increase the discovery timeout for `runTests --gtest_list_tests`
    DISCOVERY_TIMEOUT 60
)
