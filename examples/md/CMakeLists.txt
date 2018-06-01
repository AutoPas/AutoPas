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
add_test(
        NAME md-main_test
        COMMAND md-main
        CONFIGURATIONS checkExamples
)

# add the test to checkExamples
add_dependencies(checkExamples md-main)
