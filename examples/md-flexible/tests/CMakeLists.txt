if (NOT AUTOPAS_BUILD_TESTS)
    message(STATUS "AUTOPAS_BUILD_TESTS is disabled. Not building MDFlex tests!")
    return()
endif ()
message(STATUS "mdflexTests - adding tests.")

file(
    GLOB_RECURSE
    MDFlexTest_SRC
    "../src/*.cpp"
    "../src/*.h"
    "./*.cpp"
    "./*.h"
)

# remove main to avoid duplicate two main functions
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*main.cpp")
# remove files which are not needed for tests and expensive to compile
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractionsLJFunctorAVX.cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractionsLJFunctor.cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractionsLJFunctorGlobals.cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractionsLJFunctorSVE.cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractionsATFunctor.cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractionsFlopCounterFunctor.cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractionsFlopCounterFunctor3B.cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*Simulation.((cpp)|h)")

add_executable(mdFlexTests ${MDFlexTest_SRC}
        ../../../applicationLibrary/molecularDynamics/molecularDynamicsLibrary/ATLookUpTable.h
        ../../../applicationLibrary/molecularDynamics/tests/LookUpTableTests/LookUpTableTests.h
        ../../../applicationLibrary/molecularDynamics/tests/LookUpTableTests/LookUpTableTests.cpp)

target_compile_definitions(
    mdFlexTests PRIVATE
    YAMLDIRECTORY=\"${PROJECT_SOURCE_DIR}/examples/md-flexible/tests/yamlTestFiles/\"
)

include(autopas_all)

target_include_directories(
    mdFlexTests
    PUBLIC ${PROJECT_SOURCE_DIR}/tests/testAutopas ${PROJECT_SOURCE_DIR}/examples/md-flexible
)

target_link_libraries(
  mdFlexTests
  PUBLIC
    autopas
    autopasTools
    molecularDynamicsLibrary
    gmock
    yaml-cpp
    ${ALL_LIB} 
    $<$<BOOL:${MD_FLEXIBLE_USE_MPI}>:MPI::MPI_CXX>
)

if (MD_FLEXIBLE_USE_MPI)
    # use the same number of ranks as in AutoPas tests
    set_property(TARGET mdFlexTests PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${AUTOPAS_TESTS_NUM_MPI_RANKS}')
endif ()

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
gtest_discover_tests(
        mdFlexTests TEST_PREFIX "mdFlexTests/"
        PROPERTIES
            ENVIRONMENT "${LSAN_OPTIONS_STR}"
)
