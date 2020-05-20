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

add_library(autopas STATIC $<$<BOOL:${AUTOPAS_ENABLE_CUDA}>:${CU_SRC}> ${MY_SRC})

target_link_libraries(
    autopas
    PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
        spdlog::spdlog
        $<$<BOOL:${AUTOPAS_MPI}>:MPI::MPI_CXX>
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
    $<$<BOOL:${AUTOPAS_ENABLE_CUDA}>:AUTOPAS_CUDA>
    $<$<NOT:$<BOOL:${AUTOPAS_OPENMP}>>:EIGEN_DONT_PARALLELIZE>
    $<$<BOOL:${AUTOPAS_MPI}>:AUTOPAS_MPI>
)

if (AUTOPAS_ENABLE_CUDA)
    target_include_directories(autopas SYSTEM PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif ()

target_include_directories(autopas PUBLIC ../.)
