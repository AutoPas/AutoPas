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

option(AUTOPAS_LINK_WITH_IRC "Links AutoPas with irc from Intel. This is necessary on some clusters." OFF)

target_link_libraries(
    autopas
    PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
        spdlog::spdlog
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:MPI::MPI_CXX>
        $<$<BOOL:${AUTOPAS_LINK_WITH_IRC}>:irc>
        antlr4cpp
        $<$<OR:$<BOOL:${AUTOPAS_INTERNODE_TUNING}>,$<BOOL:${MD_FLEXIBLE_USE_MPI}>>:MPI::MPI_CXX>
    PRIVATE
        # harmony and Eigen3 are only needed privately when building AutoPas.
        harmony Eigen3
)

# Ompstuff needs to be here because autopas_OpenMP.cmake needs to run before this file to create to OpenMP
# target. this can be resolved by upgrading to CMake 3.13 and enforcing CMP0079.
target_compile_definitions(
    autopas
    PUBLIC
    $<$<BOOL:${AUTOPAS_OPENMP}>:AUTOPAS_OPENMP>
    $<$<NOT:$<BOOL:${AUTOPAS_OPENMP}>>:EIGEN_DONT_PARALLELIZE>
    $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:AUTOPAS_INCLUDE_MPI>
    $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:AUTOPAS_INTERNODE_TUNING>
    _USE_MATH_DEFINES
)


target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)
#target_include_directories(autopas PUBLIC "/usr/include/antlr4-runtime/")
