file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
    ../../tests/testAutopas/testingHelpers/RandomGenerator.cpp
)

add_executable(sph-traversals ${MY_SRC})

target_link_libraries(sph-traversals autopas)

# -----------------test-----------------
# add check for current target

set(ContainerList "LinkedCells" "VerletListsCells" "VerletClusterLists")

set(
    LinkedCellsTraversals
    "C01"
    "C08"
    "C18"
    "SLI"
)
set(
    VerletListsCellsTraversals
    "C01"
    ""
    "C18"
    "SLI"
)
set(
    VerletClusterListsTraversals
    ""
    ""
    ""
    ""
    "CLU"
)

set(C01NewtonModes "OFF")
set(C08NewtonModes "OFF" "ON")
set(C18NewtonModes "OFF" "ON")
set(SLINewtonModes "OFF" "ON")
set(CLUNewtonModes "OFF")

# cmake-format: off
foreach(container_type IN LISTS ContainerList)
    list(FIND ContainerList "${container_type}" containerid)
    foreach(traversal IN LISTS ${container_type}Traversals)
        list(FIND ${container_type}Traversals "${traversal}" traversalid)
        foreach(newton3 IN LISTS ${traversal}NewtonModes)
            list(FIND ${traversal}NewtonModes "${newton3}" newton3id)
            add_test(
                NAME
                    sph-traversals.test${container_type}-${traversal}-N3${newton3}
                COMMAND
                    sph-traversals 2000 10 ${containerid} 0.15 ${traversalid} ${newton3id} 0.1 5
                CONFIGURATIONS checkExamples
            )
        endforeach()
    endforeach()
endforeach()
# cmake-format: on

# add the test to checkExamples
add_dependencies(checkExamples sph-traversals)
