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
# add check for current target
add_custom_target(
        md-main_test
        COMMAND md-main > /dev/null
)
# make the test depend on current target, to ensure that the target is build
add_dependencies(md-main_test md-main)

# add the test to checkExamples
add_dependencies(checkExamples md-main_test)
