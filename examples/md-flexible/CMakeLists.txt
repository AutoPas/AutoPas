file(
        GLOB_RECURSE
        MDFlex_SRC
        "src/*.cpp"
        "src/*.h"
)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(md-flexible ${MDFlex_SRC})

target_include_directories(md-flexible PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(autopas_yaml-cpp)

option(AUTOPAS_INCLUDE_MPI "Links MPI libraries to the target. This is required if you want to use MPI." ON)

if (AUTOPAS_INTERNODE_TUNING OR AUTOPAS_INCLUDE_MPI)
    message(STATUS "AUTOPAS_INTERNODE_TUNING set, searching for mpi to link against md-flexible ...")
    find_package(MPI)

    if (NOT ${MPI_CXX_FOUND})
        message(FATAL_ERROR "cxx mpi could not be found, even though AUTOPAS_INTERNODE_TUNING was set.")
    else ()
        message(STATUS "cxx mpi found: ${MPI_CXX_COMPILER}")
        set(MD_FLEXIBLE_USE_MPI ON)
    endif ()
else ()
    message(STATUS "AUTOPAS_INTERNODE_TUNING not set, not linking MPI to md-flexible.")
endif ()

target_link_libraries(md-flexible PUBLIC autopas autopasTools yaml-cpp $<$<BOOL:${MD_FLEXIBLE_USE_MPI}>:MPI::MPI_CXX>)

# --- copy script files to build dir ---
file(
        GLOB_RECURSE SCRIPTS
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        "*.sh"
        "*.gp"
        "*.py"
        "*.yaml"
)

foreach (script ${SCRIPTS})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${script} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach (script)

# -----------------test-----------------

add_subdirectory(tests)

# add check for current target
# cmake-format: off

# randomly generated imbalanced scenario tested with LC only
# also tests the flop functor
add_test(
        NAME md-flexible.test-static
        COMMAND
        md-flexible
        --container linked
        --cutoff 1.
        --distribution-mean 5.0
        --distribution-stddeviation 2.0
        --data-layout soa
        --functor lj
        --iterations 10
        --particle-generator gauss
        --particles-total 10
        --traversal c08,sliced
        --verlet-rebuild-frequency 5
        --verlet-skin-radius 0
        --periodic false
        --deltaT 0.
        --no-end-config
        --no-progress-bar
        CONFIGURATIONS checkExamples
)

#stable, periodic particle grid tested with all configurations
add_test(
        NAME md-flexible.test-sim
        COMMAND
        md-flexible
        --container all
        --traversal all
        --cutoff 1.5
        --functor lj
        --tuning-phases 1
        --particle-generator grid
        --particles-per-dimension 10
        --particle-spacing 1.1225
        --verlet-rebuild-frequency 4
        --verlet-skin-radius 0.2
        --periodic true
        --deltaT 0.005
        --no-end-config
        --no-flops
        --no-progress-bar
        CONFIGURATIONS checkExamples
)

# The AVX example should only be tested if AVX is available.
set(HAS_AVX true)
if (NOT AUTOPAS_USE_VECTORIZATION)
    # If no vectorization is used, don't enable the example.
    message(STATUS "No vectorization used, not adding md-flexible-avx.test-unaligned to ctest.")
    set(HAS_AVX false)
elseif (AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "SSE")
    # If vectorization is set to SSE: don't enable the example
    message(STATUS "Only SSE specified, not adding md-flexible-avx.test-unaligned to ctest.")
    set(HAS_AVX false)
elseif (AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "NATIVE")
    # If Vectorization is set to native: we try to figure out the vectorization level.
    message(STATUS "Native vectorization level, trying to detect AVX capability.")
    try_compile(
            HAS_AVX
            ${CMAKE_BINARY_DIR}
            ${AUTOPAS_SOURCE_DIR}/cmake/tests/has_avx_test.cpp
            # We are abusing COMPILE_DEFINITIONS here. There is no other sane way of passing "-march=native" ...
            COMPILE_DEFINITIONS "-march=native"
            OUTPUT_VARIABLE HAS_AVX_ERROR
    )
else ()
    message(STATUS "Adding AVX example to test, as proper vectorization is manually specified.")
endif ()

if (HAS_AVX)
    message(STATUS "AVX detected. Adding AVX test.")
    #dangerous, as lc_c04_combined_SoA uses unaligned values.
    add_test(
            NAME md-flexible-avx.test-unaligned
            COMMAND
            md-flexible
            --no-end-config
            --no-flops
            --functor lennard-jones-AVX2
            --deltaT 0
            --particle-generator uniform
            --log-level debug
            --traversal lc_c04_combined_SoA
            --particles-total 71
            --no-progress-bar
            CONFIGURATIONS checkExamples
    )
else ()
    message(STATUS "AVX could not be detected. Not adding AVX test.")
endif ()

add_test(
        NAME md-flexMeasurePerf
        COMMAND measurePerf.sh md-flexible -silent
        CONFIGURATIONS checkExamples
)
# cmake-format: on

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
