FILE(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

FIND_PACKAGE(MPI)

if (NOT ${MPI_CXX_FOUND})
    message(STATUS "cxx mpi not found, not building sph-main-mpi")
    return()
else()
    message(STATUS "cxx mpi found: ${MPI_CXX_COMPILER}")
endif()

INCLUDE_DIRECTORIES(SYSTEM
        ${MPI_INCLUDE_PATH})

ADD_EXECUTABLE(sph-main-mpi
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-main-mpi
        autopas
        ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    SET_TARGET_PROPERTIES(sph-main-mpi PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    SET_TARGET_PROPERTIES(sph-main-mpi PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()


#-----------------test-----------------
# add check for current target
add_test(
        NAME sph-main-mpi.test
        COMMAND mpirun -n 4 sph-main-mpi | grep "time step" | grep -v "the time step" | tail -1 | grep -q "time step 50"
        CONFIGURATIONS checkExamples checkExamplesMPI
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples sph-main-mpi)
