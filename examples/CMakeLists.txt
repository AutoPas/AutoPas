# do nothing if nothing should be done
if (NOT AUTOPAS_BUILD_EXAMPLES)
    return()
endif ()

# md-flexible links molecularDynamicsLibrary and the sph examples link SPHLibrary. Without them the
# examples fail to compile with missing headers, so demand them here rather than at compile time.
if (NOT AUTOPAS_BUILD_MD_LIBRARY OR NOT AUTOPAS_BUILD_SPH_LIBRARY)
    message(
        FATAL_ERROR
            "AUTOPAS_BUILD_EXAMPLES=ON requires AUTOPAS_BUILD_MD_LIBRARY=ON and AUTOPAS_BUILD_SPH_LIBRARY=ON."
    )
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
