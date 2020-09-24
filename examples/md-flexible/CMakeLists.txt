file(
    GLOB_RECURSE
    MDFlex_SRC
    "src/*.cpp"
    "src/*.h"
)

add_executable(md-flexible ${MDFlex_SRC})

target_include_directories(md-flexible PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(autopas_yaml-cpp)

target_link_libraries(md-flexible PUBLIC autopas autopasTools yaml-cpp)

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

add_test(
    NAME md-flexMeasurePerf
    COMMAND measurePerf.sh md-flexible -silent
    CONFIGURATIONS checkExamples
)
# cmake-format: on

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
