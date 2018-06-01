file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(sph-main-verlet
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-main-verlet
        autopas)

#-----------------test-----------------
# add check for current target
add_custom_target(
        sph-main-verlet_test
        COMMAND sph-main-verlet | grep "time step" | tail -2 | head -1 | grep -q "time step 50"
)
# make the test depend on current target, to ensure that the target is build
add_dependencies(sph-main-verlet_test sph-main-verlet)

# add the test to checkExamples
add_dependencies(checkExamples sph-main-verlet_test)
