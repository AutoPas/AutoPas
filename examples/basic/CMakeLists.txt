file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )


ADD_EXECUTABLE(basic
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(basic
        autopas)

#-----------------test-----------------
# add check for current target
add_custom_target(
        basic_test
        COMMAND basic > /dev/null
)
# make the test depend on current target, to ensure that the target is build
add_dependencies(basic_test basic)

# add the test to checkExamples
add_dependencies(checkExamples basic_test)
