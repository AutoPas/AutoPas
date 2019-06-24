file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

find_package(MPI)

if (NOT ${MPI_CXX_FOUND})
    message(STATUS "cxx mpi not found, not building sph-main-mpi")
    return()
else ()
    message(STATUS "cxx mpi found: ${MPI_CXX_COMPILER}")
endif ()

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(sph-main-mpi ${MY_SRC})

target_link_libraries(sph-main-mpi autopas ${MPI_LIBRARIES})

if (MPI_COMPILE_FLAGS)
    set_target_properties(sph-main-mpi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(sph-main-mpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()

# -----------------test-----------------
# add check for current target
# cmake-format: off
add_test(
    NAME sph-main-mpi.test
    COMMAND mpirun -n 4 sph-main-mpi | grep "time step" | grep -v "the time step" | tail -1 | grep -q "time step 8"
    CONFIGURATIONS checkExamples checkExamplesMPI
)
# cmake-format: on

# add the executable to checkExamples as dependency
add_dependencies(checkExamples sph-main-mpi)
