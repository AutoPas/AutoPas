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
add_test(
        NAME basic.test
        COMMAND basic
        CONFIGURATIONS checkExamples
)

# add the test to checkExamples
add_dependencies(checkExamples basic)
