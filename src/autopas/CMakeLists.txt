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

target_link_libraries(
    autopas
    PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
        $<$<BOOL:${AUTOPAS_AUTO4OMP}>:auto4omp::omp>
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
        $<$<BOOL:${AUTOPAS_AUTO4OMP}>:AUTOPAS_USE_AUTO4OMP>
        $<$<BOOL:${AUTOPAS_LB4OMP}>:AUTOPAS_USE_LB4OMP>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:AUTOPAS_INCLUDE_MPI>
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:AUTOPAS_INTERNODE_TUNING>
        _USE_MATH_DEFINES
)


target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)

if (AUTOPAS_AUTO4OMP)
    #add_dependencies(autopas auto4omp)
    target_include_directories(autopas BEFORE PUBLIC ${auto4omp_SOURCE_DIR}/runtime/src) # For kmp.h, kmp_runtime.cpp
    target_include_directories(autopas BEFORE PUBLIC ${auto4omp_BINARY_DIR}/runtime/src) # For kmp_config.h
    target_link_directories(autopas BEFORE PUBLIC ${auto4omp_BINARY_DIR}/runtime/src/CMakeFiles/omp.dir) # For kmp_runtime.cpp
    target_link_directories(autopas BEFORE PUBLIC ${auto4omp_SOURCE_DIR}/runtime/src) # For kmp_runtime.cpp
    target_link_directories(autopas BEFORE PUBLIC ${auto4omp_BINARY_DIR}/runtime/src) # For kmp_runtime.cpp
    # Provide libomp.sp.5 for systems with modern OpenMP.
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp.h kmp.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_platform.h kmp_platform.h @ONLY)
    #configure_file(${auto4omp_BINARY_DIR}/runtime/src/kmp_config.h kmp_config.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_os.h kmp_os.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_safe_c_api.h kmp_safe_c_api.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_debug.h kmp_debug.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_lock.h kmp_lock.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_version.h kmp_version.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_i18n.h i18n.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_safe_c_api.h kmp_safe_c_api.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_safe_c_api.h kmp_safe_c_api.h @ONLY)
    #configure_file(${auto4omp_SOURCE_DIR}/runtime/src/kmp_safe_c_api.h kmp_safe_c_api.h @ONLY)
    add_custom_command(
            TARGET autopas PRE_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink libomp.so libomp.so.5
            WORKING_DIRECTORY "${auto4omp_BINARY_DIR}/runtime/src"
    )
endif ()
