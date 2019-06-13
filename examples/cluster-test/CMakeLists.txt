file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
    ../../tests/testAutopas/testingHelpers/RandomGenerator.cpp
)

add_executable(cluster-test ${MY_SRC})

target_link_libraries(cluster-test autopas)

# -----------------test-----------------
# add check for current target
add_test(
    NAME cluster-test.test
    COMMAND
        cluster-test
        2000 # particles
        0.12 # box size
        0.1 # skin is 0.1 * r_c
    CONFIGURATIONS checkExamples
)

# add the test to checkExamples
add_dependencies(checkExamples cluster-test)
