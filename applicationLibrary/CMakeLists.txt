# --------------------------libraries-------------------------------------------

# --- molecular dynamics ---

file(
        GLOB_RECURSE
        LIB_MD_SRC
        "molecularDynamics/molecularDynamicsLibrary/*.h"
        "molecularDynamics/molecularDynamicsLibrary/*.cpp"
)

add_library(molecularDynamicsLibrary STATIC ${LIB_MD_SRC})

# Set CMake option for md-flexible mode
set(MD_FLEXIBLE_MODE
        "SINGLESITE"
        CACHE
        STRING "Choose the mode of md-flexible to be compiled."
        )
set_property(CACHE MD_FLEXIBLE_MODE PROPERTY STRINGS "SINGLESITE;MULTISITE")

# Define numerical values for different md-flexible modes
target_compile_definitions(molecularDynamicsLibrary PUBLIC SINGLESITE=0)
target_compile_definitions(molecularDynamicsLibrary PUBLIC MULTISITE=1)

string(TOLOWER "${MD_FLEXIBLE_MODE}" MD_FLEXIBLE_MODE_lower)

if (MD_FLEXIBLE_MODE_lower MATCHES  "singlesite")
    target_compile_definitions(molecularDynamicsLibrary PUBLIC MD_FLEXIBLE_MODE=0)
elseif(MD_FLEXIBLE_MODE_lower MATCHES "multisite")
    target_compile_definitions(molecularDynamicsLibrary PUBLIC MD_FLEXIBLE_MODE=1)
endif()

target_link_libraries(molecularDynamicsLibrary PUBLIC autopas)

target_include_directories(molecularDynamicsLibrary PUBLIC ./molecularDynamics)

# --- sph ---

file(
        GLOB_RECURSE
        LIB_SPH_SRC
        "sph/SPHLibrary/*.h"
        "sph/SPHLibrary/*.cpp"
)

add_library(SPHLibrary STATIC ${LIB_SPH_SRC})

target_link_libraries(SPHLibrary PUBLIC autopas)

target_include_directories(SPHLibrary PUBLIC ./sph)

# --- proliferating cell collectives ---

file(
        GLOB_RECURSE
        LIB_PCC_SRC
        "proliferatingCellCollectives/proliferatingCellCollectivesLibrary/*.h"
        "proliferatingCellCollectives/proliferatingCellCollectivesLibrary/*.cpp"
)

add_library(proliferatingCellCollectivesLibrary STATIC ${LIB_PCC_SRC})

target_link_libraries(proliferatingCellCollectivesLibrary PUBLIC autopas molecularDynamicsLibrary)

target_include_directories(proliferatingCellCollectivesLibrary PUBLIC ./molecularDynamics/molecularDynamicsLibrary ./proliferatingCellCollectives)

# --------------------------------tests---------------------------------------

# --- molecular dynamics ---

if (AUTOPAS_BUILD_TESTS)
    message(STATUS "applicationLibrary - adding all ${MD_FLEXIBLE_MODE_lower} tests")

    if (MD_FLEXIBLE_MODE_lower MATCHES "singlesite")
        set(TEST_EXE_NAME "mdLibSingleSiteTests")
        set(TEST_EXE_UNIQUE_PATH "molecularDynamics/tests/singleSiteTests")
    elseif (MD_FLEXIBLE_MODE_lower MATCHES "multisite")
        set(TEST_EXE_NAME "mdLibMultiSiteTests")
        set(TEST_EXE_UNIQUE_PATH "molecularDynamics/tests/multiSiteTests")
    endif()

    # sanity check to prevent against globbing the whole file system
    if(NOT DEFINED TEST_EXE_NAME OR TEST_EXE_NAME STREQUAL "" OR
       NOT DEFINED TEST_EXE_UNIQUE_PATH OR TEST_EXE_UNIQUE_PATH STREQUAL "")
        message(FATAL_ERROR "Both TEST_EXE_NAME and TEST_EXE_UNIQUE_PATH must be defined and not empty.\
        This probably resulted from MD_FLEXIBLE_MODE being set wrong.\
        MD_FLEXIBLE_MODE=${MD_FLEXIBLE_MODE}")
    endif()

    file(
            GLOB_RECURSE
            TEST_SRC
            "${TEST_EXE_UNIQUE_PATH}/*.cpp"
            "${TEST_EXE_UNIQUE_PATH}/*.h"
            "molecularDynamics/tests/particlePropertiesLibraryTests/*.cpp"
            "molecularDynamics/tests/particlePropertiesLibraryTests/*.h"
            "molecularDynamics/tests/templateInstantiations/AutoPasInstantiations.cpp"
    )
    add_executable(${TEST_EXE_NAME} ${TEST_SRC} mdLibTests.cpp)

    target_include_directories(
            ${TEST_EXE_NAME}
            PUBLIC
            "${PROJECT_SOURCE_DIR}/tests/testAutopas"
            "molecularDynamics"
    )

    target_link_libraries(
            ${TEST_EXE_NAME}
            PUBLIC
            autopas autopasTools molecularDynamicsLibrary gmock
    )

    # this cmake module was only introduced in 3.10
    include(GoogleTest)
    # more robust, queries the compiled executable
    gtest_discover_tests(
            ${TEST_EXE_NAME}
            TEST_PREFIX "${TEST_EXE_NAME}/"
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
            autopas autopasTools SPHLibrary gmock
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