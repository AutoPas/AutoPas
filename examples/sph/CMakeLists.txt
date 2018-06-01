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
add_test(
        NAME sph-main_test
        COMMAND sph-main | grep "time step" | tail -2 | head -1 | grep -q "time step 50"
        CONFIGURATIONS checkExamples
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples sph-main)
