file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

add_executable(runTests ${MY_SRC})

target_compile_definitions(runTests PRIVATE)

if(AUTOPAS_ENABLE_COVERAGE)
  target_compile_options(runTests PRIVATE -O0 --coverage)

  # find required tools
  if(LCOV STREQUAL "LCOV-NOTFOUND")
    message(FATAL_ERROR "lcov not found! Please install lcov.")
  endif()
  
  find_program(GENHTML genhtml REQUIRED)
  if(GENHTML STREQUAL "GENHTML-NOTFOUND")
    message(FATAL_ERROR "genhtml not found! Please install genhtml.")
  endif()

  # add coverage target
  add_custom_target(
    coverage
    # gather data
    COMMAND
      ${LCOV} --rc lcov_branch_coverage=1 --directory . --capture --exclude */tests/* --exclude */build/* --exclude '/usr/*' --output-file coverage.info
    # generate report
    COMMAND ${GENHTML} --branch-coverage --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

target_link_libraries(
    runTests
    gcov
    autopas
    autopasTools
    molecularDynamicsLibrary
    gmock # gmock includes the gtest target, so we don't need it here.
    Eigen3 # Eigen3 is needed for the tests, but normally not required when using AutoPas.
    $<$<BOOL:${AUTOPAS_ENABLE_HARMONY}>:harmony> # ActiveHarmony is needed for the tests, but normally not required when using AutoPas.
)

include(GoogleTest)
# queries the compiled executable for tests, this requires the executable to be runnable. if you are
# cross compiling, make sure to properly set CMAKE_CROSSCOMPILING_EMULATOR.
gtest_discover_tests(
    runTests TEST_PREFIX "testAutopas/"
    # increase the discovery timeout for `runTests --gtest_list_tests`
    DISCOVERY_TIMEOUT 60
)
