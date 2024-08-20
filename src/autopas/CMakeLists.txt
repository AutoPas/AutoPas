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

option(AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        "Enables energy measurements and tuning using Intel's RAPL (Running Average Power Limit) interface for x86 processors
        running on Linux kernels."
        OFF
)

# If rules-based tuning is disabled, we must remove all files which link to antlr. As this includes generated files,
# where we can't easily change the file, we disable the files here.
if (NOT AUTOPAS_ENABLE_RULES_BASED_TUNING)
    list(FILTER MY_SRC EXCLUDE REGEX ".*parser_generated.*")
    list(FILTER MY_SRC EXCLUDE REGEX ".*RuleBasedProgramParser.*")
endif ()

add_library(autopas STATIC ${MY_SRC})

target_link_libraries(
    autopas
    PUBLIC
        $<$<NOT:$<PLATFORM_ID:Darwin>>:rt> # required for Time.h, on macOS not needed
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
        $<$<BOOL:${AUTOPAS_AUTO4OMP}>:auto4omp::omp>
        spdlog::spdlog
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:MPI::MPI_CXX>
        $<$<OR:$<BOOL:${AUTOPAS_INTERNODE_TUNING}>,$<BOOL:${MD_FLEXIBLE_USE_MPI}>>:MPI::MPI_CXX>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:MPI::MPI_CXX>
        $<$<BOOL:${AUTOPAS_ENABLE_RULES_BASED_TUNING}>:antlr4cpp>
        Eigen3
    PRIVATE
        # harmony are only needed privately when building AutoPas.
        $<$<BOOL:${AUTOPAS_ENABLE_HARMONY}>:harmony>
)

# OpenMP stuff needs to be here because autopas_OpenMP.cmake needs to run before this file to create to OpenMP
# target. this can be resolved by upgrading to CMake 3.13 and enforcing CMP0079.
target_compile_definitions(
    autopas
    PUBLIC
        $<$<BOOL:${AUTOPAS_OPENMP}>:AUTOPAS_USE_OPENMP>
        $<$<NOT:$<BOOL:${AUTOPAS_OPENMP}>>:EIGEN_DONT_PARALLELIZE>
        $<$<BOOL:${AUTOPAS_AUTO4OMP}>:AUTOPAS_USE_AUTO4OMP>
        $<$<BOOL:${AUTOPAS_LB4OMP}>:AUTOPAS_USE_LB4OMP>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:AUTOPAS_INCLUDE_MPI>
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:AUTOPAS_INTERNODE_TUNING>
        $<$<BOOL:${AUTOPAS_ENABLE_HARMONY}>:AUTOPAS_ENABLE_HARMONY>
        $<$<BOOL:${AUTOPAS_ENABLE_RULES_BASED_TUNING}>:AUTOPAS_ENABLE_RULES_BASED_TUNING>
        $<$<BOOL:${AUTOPAS_ENABLE_ENERGY_MEASUREMENTS}>:AUTOPAS_ENABLE_ENERGY_MEASUREMENTS>
        _USE_MATH_DEFINES
)


target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)

if (AUTOPAS_AUTO4OMP)
    #add_dependencies(autopas auto4omp)
    target_include_directories(autopas BEFORE PUBLIC ${auto4omp_SOURCE_DIR}/runtime/src) # For kmp.h, etc.
    target_include_directories(autopas BEFORE PUBLIC ${auto4omp_BINARY_DIR}/runtime/src) # For kmp_config.h.

    # Provide libomp.sp.5 for systems with modern OpenMP.
    add_custom_command(
            TARGET autopas PRE_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink libomp.so libomp.so.5
            WORKING_DIRECTORY "${auto4omp_BINARY_DIR}/runtime/src"
    )
endif ()
