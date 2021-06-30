file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

file(
    GLOB_RECURSE
    CU_SRC
    "*.cu"
    "*.cuh"
)

add_library(autopas STATIC ${MY_SRC})

if (AUTOPAS_INCLUDE_MPI)
    message(STATUS "AUTOPAS_INCLUDE_MPI set, searching for mpi to link against md-flexible ...")
    find_package(MPI)

    if (NOT ${MPI_CXX_FOUND})
        message(FATAL_ERROR "cxx mpi could not be found, even though AUTOPAS_INCLUDE_MPI was set.")
    else ()
        message(STATUS "cxx mpi found: ${MPI_CXX_COMPILER}")
        set(AUTOPAS_INCLUDE_MPI ON)
        add_definitions(-DAUTOPAS_INCLUDE_MPI)
    endif ()
else ()
    message(STATUS "AUTOPAS_INCLUDE_MPI not set, not linking MPI to md-flexible.")
endif ()

target_link_libraries(
    autopas
    PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
        spdlog::spdlog
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:MPI::MPI_CXX>
    PRIVATE
        # harmony and Eigen3 are only needed privately when building AutoPas.
        harmony Eigen3
)

# Ompstuff needs to be here because OpenMP.cmake needs to run before this file to create to OpenMP
# target. this can be resolved by upgrading to CMake 3.13 and enforcing CMP0079.
target_compile_definitions(
    autopas
    PUBLIC
    $<$<BOOL:${AUTOPAS_OPENMP}>:AUTOPAS_OPENMP>
    $<$<NOT:$<BOOL:${AUTOPAS_OPENMP}>>:EIGEN_DONT_PARALLELIZE>
    $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:AUTOPAS_INCLUDE_MPI>
    _USE_MATH_DEFINES
)


target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)
