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
option(AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER  "Enables fast particle buffers in all containers, thereby disabling the dynamic rebuild." OFF)
option(AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER_SELF_LIN  "Rebuild depending on amount of particles in the buffer." OFF)

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
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
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
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:AUTOPAS_INCLUDE_MPI>
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:AUTOPAS_INTERNODE_TUNING>
        $<$<BOOL:${AUTOPAS_ENABLE_HARMONY}>:AUTOPAS_ENABLE_HARMONY>
        $<$<BOOL:${AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING}>:AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING>
        $<$<BOOL:${AUTOPAS_ENABLE_ENERGY_MEASUREMENTS}>:AUTOPAS_ENABLE_ENERGY_MEASUREMENTS>
        $<$<BOOL:${AUTOPAS_ENABLE_DYNAMIC_CONTAINERS}>:AUTOPAS_ENABLE_DYNAMIC_CONTAINERS>
        $<$<BOOL:${AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER}>:AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER>
        $<$<BOOL:${AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER_SELF_LIN}>:AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER_SELF_LIN>
        _USE_MATH_DEFINES
)

target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)
