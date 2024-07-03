file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

message(STATUS "sph-main - searching for mpi...")
find_package(MPI)

if (NOT ${MPI_CXX_FOUND})
    message(STATUS "sph-main - MPI compiler NOT found. Skipping!")
    add_executable(sph-main ${MY_SRC})
    add_test (
            NAME sph-main.test
            COMMAND sph-main | grep "time step" | tail -2 | head -1 | grep -q "time step 8"
            CONFIGURATIONS checkExamples
    )
    add_dependencies(checkExamples sph-main)
else ()
    message(STATUS "sph-main - MPI compiler found: ${MPI_CXX_COMPILER}")
    add_executable(sph-main ${MY_SRC})
    target_link_libraries(sph-main autopas SPHLibrary MPI::MPI_CXX)
    # -----------------test-----------------
    # add check for current target
    # cmake-format: off
    add_test(
            NAME sph-main.test
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./sph-main-mpi | grep "time step" | grep -v "the time step" | tail -1 | grep -q "time step 8"
            CONFIGURATIONS checkExamples checkExamplesMPI
    )
    # cmake-format: on
    # add the executable to checkExamples as dependency
    add_dependencies(checkExamples sph-main)
endif ()
