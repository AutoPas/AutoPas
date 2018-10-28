# make exapmles optional. Useful when only the library is needed.
option(BUILD_EXAMPLES "Build examples." ON)

# do nothing if nothing should be done
if (NOT BUILD_EXAMPLES)
    return()
endif ()

ADD_CUSTOM_TARGET(checkExamples
        COMMAND ctest -C checkExamples -j 4)

ADD_SUBDIRECTORY(basic)
ADD_SUBDIRECTORY(md-flexible)
ADD_SUBDIRECTORY(md)
ADD_SUBDIRECTORY(sph)
ADD_SUBDIRECTORY(sph-mpi)
ADD_SUBDIRECTORY(sph-verlet)
ADD_SUBDIRECTORY(sphDiagramGeneration)
ADD_SUBDIRECTORY(sph-traversals)

