if (NOT AUTOPAS_BUILD_TESTS)
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building MDFlex tests!")
    return()
endif ()

file(
    GLOB_RECURSE
    MDFlexTest_SRC
    "../src/*.cpp"
    "../src/*.h"
    "./*.cpp"
    "./*.h"
)

# remove main
foreach (FILE_PATH ${MDFlexTest_SRC})
    string(FIND ${FILE_PATH} "main.cpp" EXCLUDE_FILE)
    if (NOT ${EXCLUDE_FILE} EQUAL -1)
        list(REMOVE_ITEM MDFlexTest_SRC ${FILE_PATH})
    endif ()
endforeach ()

add_executable(mdFlexTests ${MDFlexTest_SRC})

target_compile_definitions(
    mdFlexTests PRIVATE
    YAMLDIRECTORY=\"${PROJECT_SOURCE_DIR}/examples/md-flexible/tests/yamlTestFiles/\"
)
target_include_directories(
    mdFlexTests
    PUBLIC ${PROJECT_SOURCE_DIR}/tests/testAutopas ${PROJECT_SOURCE_DIR}/examples/md-flexible
)

target_link_libraries(
  mdFlexTests
  PUBLIC
    autopas
    autopasTools
    gmock
    yaml-cpp
    $<$<BOOL:${MD_FLEXIBLE_USE_MPI}>:MPI::MPI_CXX>
)

if (MD_FLEXIBLE_USE_MPI)
    # use the same number of ranks as in AutoPas tests
    set_property(TARGET mdFlexTests PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${AUTOPAS_TESTS_NUM_MPI_RANKS}')
endif ()

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
gtest_discover_tests(mdFlexTests TEST_PREFIX "mdFlexTests/")
