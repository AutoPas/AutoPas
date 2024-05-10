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

if (TARGET auto4omp)
    # Add the Auto4OMP dependency to AutoPas.
    add_dependencies(autopas auto4omp)

    # Include and link Auto4OMP's headers in AutoPas (e.g., omp.h).
    target_include_directories(autopas BEFORE PUBLIC "${auto4omp_BINARY_DIR}/runtime/src")
    target_link_directories(autopas BEFORE PUBLIC "${auto4omp_BINARY_DIR}/runtime/src")

    # Link Auto4OMP's built libraries.
    target_link_libraries(autopas PUBLIC auto4omp_libomp auto4omp_libomptarget)
    if (TARGET auto4omp_libomptarget_rtl_cuda)
        target_link_libraries(autopas PUBLIC auto4omp_libomptarget_rtl_cuda auto4omp_libomptarget_nvptx)
    endif ()
    if (TARGET auto4omp_libomptarget_rtl_x86_64)
        target_link_libraries(autopas PUBLIC auto4omp_libomptarget_rtl_x86_64)
    endif ()
endif ()

target_link_libraries(
    autopas
    PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
        spdlog::spdlog
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:MPI::MPI_CXX>
        antlr4cpp
        $<$<OR:$<BOOL:${AUTOPAS_INTERNODE_TUNING}>,$<BOOL:${MD_FLEXIBLE_USE_MPI}>>:MPI::MPI_CXX>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:MPI::MPI_CXX>
        Eigen3
    PRIVATE
        # harmony are only needed privately when building AutoPas.
        harmony
)

# OpenMP stuff needs to be here because autopas_OpenMP.cmake needs to run before this file to create to OpenMP
# target. this can be resolved by upgrading to CMake 3.13 and enforcing CMP0079.
target_compile_definitions(
    autopas
    PUBLIC
        $<$<BOOL:${AUTOPAS_OPENMP}>:AUTOPAS_USE_OPENMP>
        $<$<NOT:$<BOOL:${AUTOPAS_OPENMP}>>:EIGEN_DONT_PARALLELIZE>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:AUTOPAS_INCLUDE_MPI>
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:AUTOPAS_INTERNODE_TUNING>
        _USE_MATH_DEFINES
)


target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)
