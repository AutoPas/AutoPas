# make exapmles optional. Useful when only the library is needed.
option(AUTOPAS_BUILD_EXAMPLES "Build examples." ON)

# do nothing if nothing should be done
if (NOT AUTOPAS_BUILD_EXAMPLES)
    return()
endif ()

add_custom_target(
    checkExamples
    COMMAND
        ctest
        -C
        checkExamples
        -j
        4
)

add_subdirectory(md-flexible)
add_subdirectory(sph)
add_subdirectory(sph-mpi)
add_subdirectory(sphDiagramGeneration)
add_subdirectory(cluster-test)

add_subdirectory(basic-cuda)
add_subdirectory(benchmark-cuda)

add_subdirectory(fmm-integration-test)