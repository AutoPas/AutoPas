set(LCOV_ENABLED_MSG "LCOV Coverage report enabled")
set(LCOV_DISABLED_MSG "LCOV Coverage report disabled")
set(DISABLE_COVERAGE_MSG "Disabling code coverage. CMAKE_BUILD_TYPE=Debug and AUTOPAS_BUILD_TESTS=ON required.")
set(COVERAGE_REQUIREMENTS_MSG "Enable code coverage, requires GCC and CMAKE_BUILD_TYPE=Debug.")
set(NON_GNU_COMPILER_ERROR_MSG "For code coverage ensure that GCOV, LCOV and GENHTML is installed. Disabling code coverage.")

option(AUTOPAS_ENABLE_COVERAGE ${COVERAGE_REQUIREMENTS_MSG} OFF)

if (AUTOPAS_ENABLE_COVERAGE)

  find_program(GCOV_EXECUTABLE gcov)
  find_program(LCOV_EXECUTABLE lcov)
  find_program(GENHTML_EXECUTABLE genhtml)

  if (NOT GCOV_EXECUTABLE OR NOT LCOV_EXECUTABLE OR NOT GENHTML_EXECUTABLE)
    message(FATAL_ERROR ${COVERAGE_REQUIREMENTS_MSG})
  endif ()

  message(STATUS ${LCOV_ENABLED_MSG})
  if (NOT AUTOPAS_BUILD_TESTS OR NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING ${DISABLE_COVERAGE_MSG})
    set(AUTOPAS_ENABLE_COVERAGE OFF CACHE BOOL ${COVERAGE_REQUIREMENTS_MSG} FORCE)
  endif()
  if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR ${NON_GNU_COMPILER_ERROR_MSG})
    set(AUTOPAS_ENABLE_COVERAGE OFF CACHE BOOL ${COVERAGE_REQUIREMENTS_MSG} FORCE)
  endif()

  if(TARGET runTests)
    target_compile_options(runTests PRIVATE --coverage)
  endif()

  # add coverage target
  add_custom_target(
    coverage
    DEPENDS runTests
    COMMAND ${CMAKE_BINARY_DIR}/tests/testAutopas/runTests
    COMMAND ${LCOV_EXECUTABLE} --ignore-errors mismatch --rc branch_coverage=1 --directory . --capture --exclude */tests/* --exclude */build/* --exclude '/usr/*' --output-file coverage.info
    COMMAND ${GENHTML_EXECUTABLE} --branch-coverage --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
else()
  message(STATUS ${LCOV_DISABLED_MSG})
endif()