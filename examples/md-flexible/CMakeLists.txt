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
        --particles-per-dimension 10
        --particle-spacing 0.4
        --traversal c08,sliced
        --verlet-rebuild-frequency 5
        --verlet-skin-radius 0
        --periodic false
        --deltaT 0.
    CONFIGURATIONS checkExamples
)

add_test(
        NAME md-flexible.test-sim
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
        --particles-per-dimension 10
        --particle-spacing 0.4
        --traversal c08,sliced
        --verlet-rebuild-frequency 4
        --verlet-skin-radius 0.5
        --periodic true
        --deltaT 0.000000005
        CONFIGURATIONS checkExamples
)
# cmake-format: on

add_test(
    NAME md-flexMeasurePerf
    COMMAND measurePerf.sh md-flexible -silent
    CONFIGURATIONS checkExamples
)

add_test(
    NAME md-flex-VerletClusterCells-4
    COMMAND
        md-flexible
        --container
        VerletClusterCells
        --verlet-cluster-size
        4
    CONFIGURATIONS checkExamples
)

add_test(
    NAME md-flex-VerletClusterCells-32
    COMMAND
        md-flexible
        --container
        VerletClusterCells
        --verlet-cluster-size
        32
        --no-flops
    CONFIGURATIONS checkExamples
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
