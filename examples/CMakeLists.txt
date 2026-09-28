# make exapmles optional. Useful when only the library is needed.
option_conditional(AUTOPAS_BUILD_EXAMPLES
        "Build examples."
        AUTOPAS_STANDALONE_BUILD
)

# do nothing if nothing should be done
if (NOT AUTOPAS_BUILD_EXAMPLES)
    return()
endif ()

# Example tests are plain add_test() tests that need both `CONFIGURATIONS checkExamples`, which hides them from a plain
# `ctest`, and `LABELS checkExamples`, which CI selects them by.CI.
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
