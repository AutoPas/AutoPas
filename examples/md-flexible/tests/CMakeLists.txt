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

option(AUTOPAS_INCLUDE_MPI "Links MPI libraries to the target. This is required if you want to use MPI." ON)

if (AUTOPAS_INTERNODE_TUNING OR AUTOPAS_INCLUDE_MPI)
    message(STATUS "AUTOPAS_INTERNODE_TUNING set, searching for mpi to link against md-flexible ...")
    find_package(MPI)

    if (NOT ${MPI_CXX_FOUND})
        message(FATAL_ERROR "cxx mpi could not be found, even though AUTOPAS_INTERNODE_TUNING was set.")
    else ()
        message(STATUS "cxx mpi found: ${MPI_CXX_COMPILER}")
        set(MD_FLEXIBLE_USE_MPI ON)
    endif ()
else ()
    message(STATUS "AUTOPAS_INTERNODE_TUNING not set, not linking MPI to md-flexible.")
endif ()

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

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
gtest_discover_tests(mdFlexTests TEST_PREFIX "mdFlexTests/")
