file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        ../../tests/testAutopas/testingHelpers/RandomGenerator.cpp
        )

ADD_EXECUTABLE(sph-traversals
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-traversals
        autopas)


#-----------------test-----------------
# add check for current target
add_test(
        NAME sph-traversals.test-LC-C01-NoN3
        COMMAND sph-traversals 2000 10 0 0.15 0 0
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-LC-C08-NoN3
        COMMAND sph-traversals 2000 10 0 0.15 1 0
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-LC-C08-N3
        COMMAND sph-traversals 2000 10 0 0.15 1 1
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-LC-C18-NoN3
        COMMAND sph-traversals 2000 10 0 0.15 2 0
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-LC-C18-N3
        COMMAND sph-traversals 2000 10 0 0.15 2 1
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-LC-sli-NoN3
        COMMAND sph-traversals 2000 10 0 0.15 3 0
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-LC-sli-N3
        COMMAND sph-traversals 2000 10 0 0.15 3 1
        CONFIGURATIONS checkExamples
)

# verlet lists cells tests
add_test(
        NAME sph-traversals.test-VLC-C01-NoN3
        COMMAND sph-traversals 2000 10 1 0.15 0 0 0.1 5
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-VLC-C18-NoN3
        COMMAND sph-traversals 2000 10 1 0.15 2 0 0.1 5
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-VLC-C18-N3
        COMMAND sph-traversals 2000 10 1 0.15 2 1 0.1 5
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-VLC-sli-NoN3
        COMMAND sph-traversals 2000 10 1 0.15 3 0 0.1 5
        CONFIGURATIONS checkExamples
)

add_test(
        NAME sph-traversals.test-VLC-sli-N3
        COMMAND sph-traversals 2000 10 1 0.15 3 1 0.1 5
        CONFIGURATIONS checkExamples
)

# verlet cluster container
add_test(
        NAME sph-traversals.test-VerletCluster-NoN3
        COMMAND sph-traversals 2000 10 2 0.15 0 0 0.1 5
        CONFIGURATIONS checkExamples
)

# add the test to checkExamples
add_dependencies(checkExamples sph-traversals)