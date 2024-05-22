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

# Prioritize Auto4OMP's libomp by setting the LD_LIBRARY_PATH environment variable.
if (AUTOPAS_AUTO4OMP)
    if (DEFINED ENV{LD_LIBRARY_PATH} AND NOT $ENV{LD_LIBRARY_PATH} STREQUAL "")
        # If the variable's not empty, prepend Auto4OMP's libomp path.
        add_custom_command(
                TARGET autopas
                PRE_LINK COMMAND export LD_LIBRARY_PATH="${AUTOPAS_LIBOMP_PATH}:$ENV{LD_LIBRARY_PATH}"
        )
    else ()
        # Else, set it to Auto4OMP's libomp path.
        add_custom_command(
                TARGET autopas
                PRE_LINK COMMAND export LD_LIBRARY_PATH="${AUTOPAS_LIBOMP_PATH}"
        )
    endif ()
endif ()

target_link_libraries(
    autopas
    PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler

        # Standard OpenMP target, linked if (OpenMP AND NOT Auto4OMP).
        $<$<AND:$<BOOL:${AUTOPAS_OPENMP}>,$<NOT:$<BOOL:${AUTOPAS_AUTO4OMP}>>>:OpenMP::OpenMP_CXX>

        # Auto4OMP targets, linked if they exist.
        $<$<BOOL:${AUTOPAS_TARGET_omp}>:omp>
        $<$<BOOL:${AUTOPAS_TARGET_omptarget}>:omptarget>
        #$<$<BOOL:${AUTOPAS_TARGET_omptarget.rtl.cuda}>:omptarget.rtl.cuda>
        #$<$<BOOL:${AUTOPAS_TARGET_omptarget.rtl.x86_64}>:omptarget.rtl.x86_64>
        #$<$<BOOL:${AUTOPAS_TARGET_omptarget-nvptx}>:omptarget-nvptx>

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
        $<$<BOOL:${AUTOPAS_AUTO4OMP}>:AUTOPAS_LIBOMP_PATH="${AUTOPAS_LIBOMP_PATH}">
        _USE_MATH_DEFINES
)


target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)
