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
add_test(
        NAME sph-diagram-generation.test_direct_sum_density
        COMMAND sph-diagram-generation 1000 2 1 0
        CONFIGURATIONS checkExamples checkExamplesSPHDia
)
add_test(
        NAME sph-diagram-generation.test_linked_cells_density
        COMMAND sph-diagram-generation 1000 2 0 0
        CONFIGURATIONS checkExamples checkExamplesSPHDia
)
add_test(
        NAME sph-diagram-generation.test_direct_sum_hydro
        COMMAND sph-diagram-generation 1000 2 1 1
        CONFIGURATIONS checkExamples checkExamplesSPHDia
)
add_test(
        NAME sph-diagram-generation.test_linked_cells_hydro
        COMMAND sph-diagram-generation 1000 2 0 1
        CONFIGURATIONS checkExamples checkExamplesSPHDia
)


# add the executable to checkExamples as dependency
add_dependencies(checkExamples sph-diagram-generation)

#sph-diagram-generation ${Mols[$i]} ${Reps[$i]} ${iCont} 0