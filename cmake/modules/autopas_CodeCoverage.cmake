set(LCOV_ENABLED_MSG "LCOV Coverage report enabled")
set(LCOV_DISABLED_MSG "LCOV Coverage report disabled")
set(LCOV_REQUIRED_VERSION_MSG "Disabling code coverage. LCOV version 2.0 or higher is required.")
set(DISABLE_COVERAGE_MSG "Disabling code coverage. CMAKE_BUILD_TYPE=Debug and AUTOPAS_BUILD_TESTS=ON required.")
set(COVERAGE_REQUIREMENTS_MSG "Enable code coverage, requires GCC and CMAKE_BUILD_TYPE=Debug.")
set(NON_GNU_COMPILER_ERROR_MSG "For code coverage ensure that GCOV, LCOV and GENHTML is installed. Disabling code coverage.")

option(AUTOPAS_ENABLE_COVERAGE ${COVERAGE_REQUIREMENTS_MSG} OFF)

if (AUTOPAS_ENABLE_COVERAGE)

  # search for required tools (lcov >= 2.1 is required)
  find_program(GCOV_EXECUTABLE gcov)
  find_program(LCOV_EXECUTABLE lcov)
  find_program(GENHTML_EXECUTABLE genhtml)

  if (NOT GCOV_EXECUTABLE OR NOT LCOV_EXECUTABLE OR NOT GENHTML_EXECUTABLE)
    message(FATAL_ERROR ${COVERAGE_REQUIREMENTS_MSG})
    else()
    message(STATUS "GCOV found at: "${GCOV_EXECUTABLE})
    message(STATUS "LCOV found at: "${LCOV_EXECUTABLE})
    message(STATUS "GENHATML found at: "${GENHTML_EXECUTABLE})
  endif ()

  # check for lcov >= 2.0
  # Execute 'lcov --version' and capture the output
  execute_process(
      COMMAND ${LCOV_EXECUTABLE} --version
      OUTPUT_VARIABLE LCOV_VERSION_OUTPUT
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  
  # Check if the output contains "v2." indicating version 2 or higher
  if(NOT LCOV_VERSION_OUTPUT MATCHES "LCOV version ([2-9]|[1-9][0-9]+)\\.([0-9]+)")
    message(FATAL_ERROR ${LCOV_REQUIRED_VERSION_MSG})
  endif()

  message(STATUS ${LCOV_ENABLED_MSG})
  
  # Debug mode has to be enabled for creating coverage reports
  if (NOT AUTOPAS_BUILD_TESTS OR NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(FATAL_ERROR ${DISABLE_COVERAGE_MSG})
    set(AUTOPAS_ENABLE_COVERAGE OFF CACHE BOOL ${COVERAGE_REQUIREMENTS_MSG} FORCE)
  endif()

  # GCOV comes with GCC and only works with GCC as compiler
  if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR ${NON_GNU_COMPILER_ERROR_MSG})
    set(AUTOPAS_ENABLE_COVERAGE OFF CACHE BOOL ${COVERAGE_REQUIREMENTS_MSG} FORCE)
  endif()

  # set the coverage flag for the test executable
  if(TARGET runTests)
    target_compile_options(runTests PRIVATE --coverage)
  endif()

  # add coverage target
  add_custom_target(
    coverage
    # test have to be built before
    DEPENDS runTests
    # first run the tests
    COMMAND ${CMAKE_BINARY_DIR}/tests/testAutopas/runTests
    # enable branch coverage and exclude uninteresting paths from coverage reporting.
    # If we do not use --ignore-errors mismatch the coverage report aborts, since there seem to be some inconsistent entries.
    COMMAND ${LCOV_EXECUTABLE} --ignore-errors mismatch --ignore-errors inconsistent --rc branch_coverage=1 --directory . --capture --exclude '*/tests/*' --exclude '/usr/*' --exclude '*/build/*' --output-file coverage.info
    # Convert .info output to html report
    COMMAND ${GENHTML_EXECUTABLE} --branch-coverage --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
else()
  message(STATUS ${LCOV_DISABLED_MSG})
endif()