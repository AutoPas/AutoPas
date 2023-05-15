file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

add_executable(sph-diagram-generation ${MY_SRC})

target_link_libraries(sph-diagram-generation autopas autopasTools SPHLibrary)

# -----------------test-----------------
# add check for current target
# cmake-format: off
add_test(
    NAME sph-diagram-generation.test_direct_sum_density
    COMMAND sph-diagram-generation 1000 2 directSum 0
    CONFIGURATIONS checkExamples checkExamplesSPHDia
)
add_test(
    NAME sph-diagram-generation.test_linked_cells_density
    COMMAND sph-diagram-generation 1000 2 linkedCells 0
    CONFIGURATIONS checkExamples checkExamplesSPHDia
)
add_test(
    NAME sph-diagram-generation.test_direct_sum_hydro
    COMMAND sph-diagram-generation 1000 2 directSum 1
    CONFIGURATIONS checkExamples checkExamplesSPHDia
)
add_test(
    NAME sph-diagram-generation.test_linked_cells_hydro
    COMMAND sph-diagram-generation 1000 2 linkedCells 1
    CONFIGURATIONS checkExamples checkExamplesSPHDia
)
# cmake-format: on

# add the executable to checkExamples as dependency
add_dependencies(checkExamples sph-diagram-generation)

# sph-diagram-generation ${Mols[$i]} ${Reps[$i]} ${iCont} 0
