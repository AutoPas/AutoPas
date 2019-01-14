file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        ../../tests/testAutopas/testingHelpers/RandomGenerator.cpp
        )

ADD_EXECUTABLE(cluster-test
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(cluster-test
        autopas)

#-----------------test-----------------
# add check for current target

add_test(
        NAME cluster-test.test
        # 2000 particles, box size is 0.12, skin is 0.1 * r_c
        COMMAND cluster-test 2000 0.12 0.1
        CONFIGURATIONS checkExamples
)

# add the test to checkExamples
add_dependencies(checkExamples cluster-test)
