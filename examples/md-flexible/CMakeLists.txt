file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
    "../../tests/testAutopas/testingHelpers/RandomGenerator.cpp"
)

add_executable(md-flexible ${MY_SRC})

include_directories(${PROJECT_SOURCE_DIR}/libs/config4cpp/include)
link_directories(${PROJECT_SOURCE_DIR}/libs/config4cpp/lib)
target_link_libraries(md-flexible autopas config4cpp)

# --- copy script files to build dir ---
file(GLOB_RECURSE SCRIPTS RELATIVE ${CMAKE_SOURCE_DIR} "*.sh" "*.gp" "*.py")

foreach (script ${SCRIPTS})
    configure_file(${CMAKE_SOURCE_DIR}/${script} ${CMAKE_BINARY_DIR}/${script} COPYONLY)
endforeach (script)

# -----------------test-----------------
# add check for current target
# cmake-format: off
add_test(
    NAME md-flexible.test
    COMMAND
        md-flexible
        --container linked
        --cutoff 1.
        --distribution-mean 5.0
        --distribution-stddeviation 2.0
        --data-layout soa
        --functor lj
        --iterations 10
        --particles-generator gauss
        --particles-per-dimension 10
        --particle-spacing 0.4
        --traversal c08,sliced
        --verlet-rebuild-frequency 5
        --verlet-skin-radius 0
    CONFIGURATIONS checkExamples
)
# cmake-format: on

add_test(
    NAME md-flexMeasurePerf
    COMMAND ${PROJECT_SOURCE_DIR}/examples/md-flexible/measurePerf.sh md-flexible -silent
    CONFIGURATIONS checkExamples
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
