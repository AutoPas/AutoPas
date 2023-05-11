file(
        GLOB_RECURSE
        LIB_SRC
        "molecularDynamics/molecularDynamicsLibrary/*.h"
)

add_library(molecularDynamicsLibrary STATIC ${LIB_SRC})

target_link_libraries(molecularDynamicsLibrary PUBLIC autopas)

target_include_directories(molecularDynamicsLibrary PUBLIC ./molecularDynamics)

# --------------------------------tests---------------------------------------

if (AUTOPAS_BUILD_TESTS)
    file(
            GLOB_RECURSE
            SS_TEST_SRC
            "molecularDynamics/singleSiteTests/*.cpp"
            "molecularDynamics/singleSiteTests/*.h"
    )

    add_executable(mdLibSingleSiteTests ${SS_TEST_SRC})

    target_compile_definitions(
            mdLibSingleSiteTests PRIVATE
            MD_FLEXIBLE_MODE=0 # 0 corresponds to single-site
    )

    target_include_directories(
            mdLibSingleSiteTests PUBLIC
            ${PROJECT_SOURCE_DIR}/tests/testAutopas molecularDynamics
    )

    target_link_libraries(
            mdLibSingleSiteTests PUBLIC
            autopas autopasTools gmock
    )

    file(
            GLOB_RECURSE
            MS_TEST_SRC
            "molecularDynamics/multiSiteTests/*.cpp"
            "molecularDynamics/multiSiteTests/*.h"
    )

    add_executable(mdLibMultiSiteTests ${MS_TEST_SRC})

    target_compile_definitions(
            mdLibMultiSiteTests PRIVATE
            MD_FLEXIBLE_MODE=1 # 1 corresponds to multi-site
    )

    target_include_directories(
            mdLibMultiSiteTests PUBLIC
            ${PROJECT_SOURCE_DIR}/tests/testAutopas molecularDynamics
    )

    target_link_libraries(
            mdLibMultiSiteTests PUBLIC
            autopas autopasTools gmock
    )

    # this cmake module was only introduced in 3.10
    include(GoogleTest)
    # more robust, queries the compiled executable
    gtest_discover_tests(
            mdLibSingleSiteTests TEST_PREFIX "mdLibSingleSiteTests/"
            PROPERTIES
            ENVIRONMENT "${LSAN_OPTIONS_STR}"
    )
    gtest_discover_tests(
            mdLibMultiSiteTests TEST_PREFIX "mdLibMultiSiteTests/"
            PROPERTIES
            ENVIRONMENT "${LSAN_OPTIONS_STR}"
    )

else()
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building Molecular Dynamics Library tests!")
endif()