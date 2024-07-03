file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)
option(SPH_USE_MPI "Enable MPI parallelism for sph-main." OFF)

message(STATUS "sph-main - searching for mpi...")
find_package(MPI)

if (NOT ${MPI_CXX_FOUND})
    message(STATUS "sph-main - MPI compiler NOT found. Skipping!")
    return()
else ()
    message(STATUS "sph-main - MPI compiler found: ${MPI_CXX_COMPILER}")
endif ()

add_executable(sph-main ${MY_SRC})
target_link_libraries(sph-main autopas SPHLibrary $<$<BOOL:${SPH_USE_MPI}>:MPI::MPI_CXX>)

if ($<$<BOOL:${SPH_USE_MPI})
    #-----------------test-----------------
    # add check for current target
    # cmake-format: off
    add_test(
            NAME sph-main.test
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./sph-main | grep "time step" | grep -v "the time step" | tail -1 | grep -q "time step 8"
            CONFIGURATIONS checkExamples checkExamplesMPI
    )
    # cmake-format: on
    # add the executable to checkExamples as dependency
    add_dependencies(checkExamples sph-main)
else()
    add_test (
            NAME sph-main.test
            COMMAND sph-main | grep "time step" | tail -2 | head -1 | grep -q "time step 8"
            CONFIGURATIONS checkExamples
    )
    add_dependencies(checkExamples sph-main)
endif()

