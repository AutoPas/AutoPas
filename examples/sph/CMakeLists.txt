file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(sph-main
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-main
        autopas)


#-----------------test-----------------
# add check for current target
add_custom_target(
        sph-main_test
        COMMAND sph-main | grep "time step" | tail -2 | head -1 | grep -q "time step 50"
)
# make the test depend on current target, to ensure that the target is build
add_dependencies(sph-main_test sph-main)

# add the test to checkExamples
add_dependencies(checkExamples sph-main_test)
