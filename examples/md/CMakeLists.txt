file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(md-main
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(md-main
        autopas)

#-----------------test-----------------
# default execttion
add_test(
        NAME md-main.test-default
        COMMAND md-main
        CONFIGURATIONS checkExamples
)
# linked cells
add_test(
        NAME md-main.test-linked-cells
        COMMAND md-main 0 1000 2
        CONFIGURATIONS checkExamples
)
# linked cells soa
add_test(
        NAME md-main.test-linked-cells-soa
        COMMAND md-main 0 1000 2 soa
        CONFIGURATIONS checkExamples
)
# verlet
add_test(
        NAME md-main.test-verlet
        COMMAND md-main 2 1000 4 2 0.1
        CONFIGURATIONS checkExamples
)
# verlet soa
add_test(
        NAME md-main.test-verlet-soa
        COMMAND md-main 2 1000 4 2 0.1 soa
        CONFIGURATIONS checkExamples
)
# direct sum
add_test(
        NAME md-main.test-direct-sum
        COMMAND md-main 1 1000 2
        CONFIGURATIONS checkExamples
)


# add the test to checkExamples
add_dependencies(checkExamples md-main)
