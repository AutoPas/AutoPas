file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(sph-diagram-generation
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-diagram-generation
        autopas)


#-----------------test-----------------
# add check for current target
add_custom_target(
        sph-diagram-generation_test
        COMMAND sph-diagram-generation 1000 2 0 0 > /dev/null
        COMMAND sph-diagram-generation 1000 2 1 0 > /dev/null
)
# make the test depend on current target, to ensure that the target is build
add_dependencies(sph-diagram-generation_test sph-diagram-generation)

# add the test to checkExamples
add_dependencies(checkExamples sph-diagram-generation_test)

#sph-diagram-generation ${Mols[$i]} ${Reps[$i]} ${iCont} 0