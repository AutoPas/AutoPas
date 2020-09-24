file(
        GLOB_RECURSE
        MDFlex_SRC
        "src/*.cpp"
        "src/*.h"
)

message(STATUS "searching for mpi...")
find_package(MPI)

if (NOT ${MPI_CXX_FOUND})
    message(STATUS "cxx mpi not found, not building sph-main-mpi")
    return()
else ()
    message(STATUS "cxx mpi found: ${MPI_CXX_COMPILER}")
endif ()

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(md-flexible ${MDFlex_SRC})

target_include_directories(md-flexible PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(autopas_yaml-cpp)

target_link_libraries(md-flexible PUBLIC autopas autopasTools yaml-cpp MPI::MPI_CXX)

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
        CONFIGURATIONS checkExamples
)

set(HAS_AVX true)
if (NOT AUTOPAS_USE_VECTORIZATION)
    message(STATUS "No vectorization used, not adding md-flexible-avx.test-unaligned to ctest.")
    set(HAS_AVX false)
elseif(AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "SSE")
    message(STATUS "Only SSE specified, not adding md-flexible-avx.test-unaligned to ctest.")
    set(HAS_AVX false)
elseif(AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "NATIVE")
    message(STATUS "Native vectorization level, trying to detect AVX capability.")
    if(UNIX)
        message(STATUS "Unix")
        EXEC_PROGRAM(cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO)
        STRING(REGEX REPLACE "^.*(avx).*$" "\\1" AVX_THERE ${CPUINFO})
        STRING(COMPARE EQUAL "avx" "${AVX_THERE}" HAS_AVX)
        message(STATUS "detected avx on Unix? ${HAS_AVX}")
    else()
        message(STATUS "no Unix")
        # set to false, as we don't know...
        set(HAS_AVX false)
    endif()
else()
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
            CONFIGURATIONS checkExamples
    )
endif ()

add_test(
        NAME md-flexMeasurePerf
        COMMAND measurePerf.sh md-flexible -silent
        CONFIGURATIONS checkExamples
)
# cmake-format: on

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
