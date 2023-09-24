option(AUTOPAS_CODE_COVERAGE "Code Coverage" OFF)

if (AUTOPAS_CODE_COVERAGE MATCHES ON)
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        message(STATUS "ENABLING CODE COVERAGE TARGET")
        include(CodeCoverage)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_COVERAGE} ${CMAKE_CXX_FLAGS_DEBUG}")
        setup_target_for_coverage(${PROJECT_NAME}_coverage runTests coverage)
        setup_target_for_coverage_cobertura(
            ${PROJECT_NAME}_cobertura
            runTests
            coverage
            --gtest_output=xml:coverage.junit.xml
            .*gtest.*
        )
    else ()
        message(FATAL_ERROR "Code Coverage only works with CMAKE_BUILD_TYPE = Debug")
    endif ()
endif ()
