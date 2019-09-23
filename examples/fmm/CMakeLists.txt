file(
        GLOB_RECURSE
        MY_SRC
        "*.cpp"
        "*.h"
        "../../tests/testAutopas/testingHelpers/RandomGenerator.cpp"
)

add_executable(fmm ${MY_SRC})

target_link_libraries(fmm autopas)

# add the test to checkExamples
add_dependencies(checkExamples fmm)