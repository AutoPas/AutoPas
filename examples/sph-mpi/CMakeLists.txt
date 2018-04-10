FILE(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

FIND_PACKAGE(MPI)

if (NOT ${MPI_CXX_FOUND})
    message(AUTHOR_WARNING "cxx mpi not found")
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