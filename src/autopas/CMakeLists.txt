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


set(AUTOPAS_ENABLE_DYNAMIC_CONTAINERS "ON"
    CACHE STRING "Dynamic container rebuilding: OFF, ON, or FPB")
set_property(CACHE AUTOPAS_ENABLE_DYNAMIC_CONTAINERS PROPERTY STRINGS OFF ON FPB)

if(NOT AUTOPAS_ENABLE_DYNAMIC_CONTAINERS MATCHES "^(OFF|ON|FPB)$")
  message(FATAL_ERROR
    "AUTOPAS_ENABLE_DYNAMIC_CONTAINERS must be OFF, ON, or FPB "
    "(got '${AUTOPAS_ENABLE_DYNAMIC_CONTAINERS}')")
endif()

# NOTE: AUTOPAS_ENABLE_DYNAMIC_CONTAINERS must be defined for every non-OFF mode;
# shared code guards on it. Mode-specific code must guard on exactly one of
# AUTOPAS_DYNAMIC_CONTAINERS_{STANDARD,FPB}.
set(_dyn_defs AUTOPAS_ENABLE_DYNAMIC_CONTAINERS)   # ← base macro, added for BOTH ON and FPB
if(AUTOPAS_ENABLE_DYNAMIC_CONTAINERS STREQUAL "ON")
  list(APPEND _dyn_defs AUTOPAS_DYNAMIC_CONTAINERS_STANDARD)
elseif(AUTOPAS_ENABLE_DYNAMIC_CONTAINERS STREQUAL "FPB")
  list(APPEND _dyn_defs AUTOPAS_DYNAMIC_CONTAINERS_FPB)
endif()

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
        hwy
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
        ${_dyn_defs}
        _USE_MATH_DEFINES
)

target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)
