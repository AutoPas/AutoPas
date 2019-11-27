file(
    GLOB_RECURSE
    MDFlex_SRC
    "*.cpp"
    "*.h"
    "${PROJECT_SOURCE_DIR}/tests/testAutopas/testingHelpers/RandomGenerator.cpp"
)

# remove all files from the md-flexible/tests dir
foreach (FILE_PATH ${MDFlex_SRC})
    string(FIND ${FILE_PATH} "md-flexible/tests/" EXCLUDE_FILE)
    if (NOT ${EXCLUDE_FILE} EQUAL -1)
        list(REMOVE_ITEM MDFlex_SRC ${FILE_PATH})
    endif ()
endforeach ()

add_executable(md-flexible ${MDFlex_SRC})

target_include_directories(md-flexible PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(autopas_yaml-cpp)

target_link_libraries(md-flexible PUBLIC autopas yaml-cpp)

# --- copy script files to build dir ---
file(
    GLOB_RECURSE SCRIPTS
    RELATIVE ${CMAKE_SOURCE_DIR}
    "*.sh"
    "*.gp"
    "*.py"
    "*.yaml"
)

foreach (script ${SCRIPTS})
    configure_file(${CMAKE_SOURCE_DIR}/${script} ${CMAKE_BINARY_DIR}/examples/md-flexible COPYONLY)
endforeach (script)

# -----------------test-----------------

add_subdirectory(tests)

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
        --particle-generator gauss
        --particles-per-dimension 10
        --particle-spacing 0.4
        --traversal c08,sliced
        --verlet-rebuild-frequency 5
        --verlet-skin-radius 0
        --periodic false
    CONFIGURATIONS checkExamples
)
# cmake-format: on

add_test(
    NAME md-flexMeasurePerf
    COMMAND measurePerf.sh md-flexible -silent
    CONFIGURATIONS checkExamples
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
