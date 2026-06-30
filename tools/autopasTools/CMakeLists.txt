file(GLOB MY_ROOT_HEADERS "*.h")
file(GLOB_RECURSE MY_GENERATOR_SRC "generators/*.cpp" "generators/*.h")
set(MY_SRC ${MY_ROOT_HEADERS} ${MY_GENERATOR_SRC})

add_library(autopasTools STATIC ${MY_SRC})

target_link_libraries(autopasTools PUBLIC autopas)

target_include_directories(autopasTools PUBLIC ..)

# --------------------------------tests---------------------------------------

if (AUTOPAS_BUILD_TESTS)
    file(
            GLOB_RECURSE
            AUTOPAS_TOOLS_TEST_SRC
            "tests/*.cpp"
            "tests/*.h"
    )

    add_executable(autopasToolsTests ${AUTOPAS_TOOLS_TEST_SRC} autopasToolsTests.cpp)

    target_include_directories(
            autopasToolsTests
            PUBLIC
            "${PROJECT_SOURCE_DIR}/tests/testAutopas"
            "tests"
    )

    target_link_libraries(
            autopasToolsTests
            PUBLIC
            autopas autopasTools molecularDynamicsLibrary gmock
    )

    include(GoogleTest)
    gtest_discover_tests(
            autopasToolsTests
            TEST_PREFIX "autopasToolsTests/"
            DISCOVERY_MODE PRE_TEST
            PROPERTIES
            ENVIRONMENT "${LSAN_OPTIONS_STR}"
    )
else()
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building autopasTools tests!")
endif()
