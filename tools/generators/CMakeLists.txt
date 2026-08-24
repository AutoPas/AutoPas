
file(GLOB_RECURSE MY_GENERATOR_SRC "src/*.cpp" "src/*.h")

add_library(generators STATIC ${MY_GENERATOR_SRC})

target_link_libraries(generators PUBLIC autopas)

target_include_directories(generators PUBLIC ..)

# --------------------------------tests---------------------------------------

if (AUTOPAS_BUILD_TESTS)
    file(
            GLOB_RECURSE
            GENERATORS_TOOLS_TEST_SRC
            "tests/*.cpp"
            "tests/*.h"
    )

    add_executable(generatorsTests ${GENERATORS_TOOLS_TEST_SRC} generatorsTests.cpp)

    target_include_directories(
            generatorsTests
            PUBLIC
            "${PROJECT_SOURCE_DIR}/tests/testAutopas"
            "tests"
    )

    target_link_libraries(
            generatorsTests
            PUBLIC
            autopas generators molecularDynamicsLibrary gmock
    )

    include(GoogleTest)
    gtest_discover_tests(
            generatorsTests
            TEST_PREFIX "generatorsTests/"
            DISCOVERY_MODE PRE_TEST
            PROPERTIES
            ENVIRONMENT "${LSAN_OPTIONS_STR}"
    )
else ()
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building generators tests!")
endif ()
