# --------------------------libraries-------------------------------------------

# --- molecular dynamics ---

file(
        GLOB_RECURSE
        LIB_MD_SRC
        "molecularDynamics/molecularDynamicsLibrary/*.h"
)

add_library(molecularDynamicsLibrary STATIC ${LIB_MD_SRC})

target_link_libraries(molecularDynamicsLibrary PUBLIC autopas)

target_include_directories(molecularDynamicsLibrary PUBLIC ./molecularDynamics)

# --- sph ---

file(
        GLOB_RECURSE
        LIB_SPH_SRC
        "sph/SPHLibrary/*.h"
)

add_library(SPHLibrary STATIC ${LIB_SPH_SRC})

target_link_libraries(SPHLibrary PUBLIC autopas)

target_include_directories(SPHLibrary PUBLIC ./sph)

# --------------------------------tests---------------------------------------

# --- molecular dynamics ---

if (AUTOPAS_BUILD_TESTS)
    message(STATUS "applicationLibrary - adding all tests")
    file(
            GLOB_RECURSE
            SS_TEST_SRC
            "molecularDynamics/singleSiteTests/*.cpp"
            "molecularDynamics/singleSiteTests/*.h"
            "molecularDynamics/particlePropertiesLibraryTests/*.cpp"
            "molecularDynamics/particlePropertiesLibraryTests/*.h"
            "molecularDynamics/templateInstatiations/AutoPasInstantiations.cpp"
    )

    add_executable(mdLibSingleSiteTests ${SS_TEST_SRC} mdLibTests.cpp)

    target_compile_definitions(
            mdLibSingleSiteTests PRIVATE
            MD_FLEXIBLE_MODE=0 # set MD_FLEXIBLE_MODE to single-site
            # Tests need definitions of SINGLESITE and MULTISITE so it can correctly compile single-site tests
            SINGLESITE=0
            MULTISITE=1
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
            "molecularDynamics/particlePropertiesLibraryTests/*.cpp"
            "molecularDynamics/particlePropertiesLibraryTests/*.h"
            "molecularDynamics/templateInstatiations/AutoPasInstantiations.cpp"
    )

    add_executable(mdLibMultiSiteTests ${MS_TEST_SRC} mdLibTests.cpp)

    target_compile_definitions(
            mdLibMultiSiteTests PRIVATE
            MD_FLEXIBLE_MODE=1  # set MD_FLEXIBLE_MODE to multi-site
            # Tests need definitions of SINGLESITE and MULTISITE so it can correctly compile multi-site tests
            SINGLESITE=0
            MULTISITE=1
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

# --- SPH ---
    file(
            GLOB_RECURSE
            SPH_TEST_SRC
            "sph/tests/*.cpp"
            "sph/tests/*.h"
    )

    add_executable(sphLibTests ${SPH_TEST_SRC} sphLibTests.cpp)

    target_include_directories(
            sphLibTests PUBLIC
            ${PROJECT_SOURCE_DIR}/tests/testAutopas sph
    )

    target_link_libraries(
            sphLibTests PUBLIC
            autopas autopasTools gmock
    )

    # this cmake module was only introduced in 3.10
    include(GoogleTest)
    # more robust, queries the compiled executable
    gtest_discover_tests(
            sphLibTests TEST_PREFIX "sphLibTests/"
            PROPERTIES
            ENVIRONMENT "${LSAN_OPTIONS_STR}"
    )

else()
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building Application Libraries tests!")
endif()