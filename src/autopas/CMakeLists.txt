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

option(AUTOPAS_ENABLE_DYNAMIC_CONTAINERS  "Enables dynamic rebuilding of all containers." ON)

# If rules-based tuning and fuzzy tuning are disabled, we must remove all files which link to antlr. As this includes generated files,
# where we can't easily change the file, we disable the files here.
if (NOT AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING)
    list(FILTER MY_SRC EXCLUDE REGEX ".*parser_generated.*")
    list(FILTER MY_SRC EXCLUDE REGEX ".*RuleBasedProgramParser.*")
endif ()

add_library(autopas STATIC ${MY_SRC})

target_link_libraries(
    autopas
    PUBLIC
        $<$<NOT:$<PLATFORM_ID:Darwin>>:rt> # required for Time.h, on macOS not needed
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_AUTO4OMP}>:auto4omp::omp>
        spdlog::spdlog
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:MPI::MPI_CXX>
        $<$<OR:$<BOOL:${AUTOPAS_INTERNODE_TUNING}>,$<BOOL:${MD_FLEXIBLE_USE_MPI}>>:MPI::MPI_CXX>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:MPI::MPI_CXX>
        $<$<BOOL:${AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING}>:antlr4cpp>
        $<$<BOOL:${AUTOPAS_ENABLE_ENERGY_MEASUREMENTS}>:pmt>
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
        $<$<BOOL:${AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING}>:AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING>
        $<$<BOOL:${AUTOPAS_ENABLE_ENERGY_MEASUREMENTS}>:AUTOPAS_ENABLE_ENERGY_MEASUREMENTS>
        $<$<BOOL:${AUTOPAS_ENABLE_DYNAMIC_CONTAINERS}>:AUTOPAS_ENABLE_DYNAMIC_CONTAINERS>
        _USE_MATH_DEFINES
)

# -fopenmp enables pragma recognition and __kmpc_* code generation (required for LB4OMP dispatch).
# Must be in compile_options (not link_options) so the compiler emits __kmpc_dispatch_* calls without
# also implicitly linking system libomp — only LB4OMP's libomp.so (from auto4omp::omp) is linked.
# -fopenmp-version=45 restricts to OMP 4.5 to avoid calling functions absent from LB4OMP's older LLVM base.
target_compile_options(autopas PUBLIC
        $<$<BOOL:${AUTOPAS_OPENMP}>:-fopenmp -fopenmp-version=45>)

target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)

if (AUTOPAS_AUTO4OMP)
    target_include_directories(autopas BEFORE PUBLIC ${auto4omp_SOURCE_DIR}/runtime/src)
    target_include_directories(autopas BEFORE PUBLIC ${auto4omp_BINARY_DIR}/runtime/src)

    # Ensure Auto4OMP's libomp.so is searched before the system's at both link and runtime
    target_link_directories(autopas BEFORE PUBLIC ${auto4omp_BINARY_DIR}/runtime/src)
    set_target_properties(autopas PROPERTIES BUILD_RPATH "${auto4omp_BINARY_DIR}/runtime/src")
endif ()

