file(
        GLOB_RECURSE
        MY_SRC
        "*.cpp"
        "*.h"
        "../../tests/testAutopas/testingHelpers/RandomGenerator.cpp"
)

add_executable(fmm-integration-test ${MY_SRC} main.cpp)

target_link_libraries(fmm-integration-test autopas)


# add the executable to checkExamples as dependency
add_dependencies(checkExamples fmm-integration-test)
