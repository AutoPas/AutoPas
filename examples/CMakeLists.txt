# make exapmles optional. Useful when only the library is needed.
option_conditional(AUTOPAS_BUILD_EXAMPLES
        "Build examples."
        AUTOPAS_STANDALONE_BUILD
)

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
add_subdirectory(md-julia)
# add_subdirectory(sph)
# add_subdirectory(sph-mpi)
# add_subdirectory(sphDiagramGeneration)
