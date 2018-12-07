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
add_test(
        NAME sph-main-verlet.test
        COMMAND sph-main-verlet | grep "time step" | tail -2 | head -1 | grep -q "time step 8"
        CONFIGURATIONS checkExamples
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples sph-main-verlet)
