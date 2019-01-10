file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        ../../tests/testAutopas/testingHelpers/RandomGenerator.cpp
        )

ADD_EXECUTABLE(sph-traversals
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-traversals
        autopas)


#-----------------test-----------------
# add check for current target


set(ContainerList "LinkedCells" "VerletListsCells" "VerletClusterLists")

set(LinkedCellsTraversals        "C01" "C08" "C18" "SLI")
set(VerletListsCellsTraversals   "C01" ""    "C18" "SLI")
set(VerletClusterListsTraversals "C01" )

set(C01NewtonModes "OFF")
set(C08NewtonModes "OFF" "ON")
set(C18NewtonModes "OFF" "ON")
set(SLINewtonModes "OFF" "ON")

foreach(container_type IN LISTS ContainerList)
    list(FIND ContainerList "${container_type}" containerid)
    foreach(traversal IN LISTS ${container_type}Traversals)
        list(FIND ${container_type}Traversals "${traversal}" traversalid)
        foreach(newton3 IN LISTS ${traversal}NewtonModes)
            list(FIND ${traversal}NewtonModes "${newton3}" newton3id)
            message("${container_type}=${containerid} ${traversal}=${traversalid} ${newton3}=${newton3id}")
            add_test(
                    NAME sph-traversals.test${container_type}-${traversal}-N3${newton3}
                    COMMAND sph-traversals 2000 10 ${containerid} 0.15 ${traversalid} ${newton3id} 0.1 5
                    CONFIGURATIONS checkExamples
            )
        endforeach()
    endforeach()
endforeach()


# add the test to checkExamples
add_dependencies(checkExamples sph-traversals)