file(
        GLOB_RECURSE
        MY_SRC
        "*.cpp"
        "*.h"
        "../../tools/autopasTools/generators/RandomGenerator.cpp"
)

add_executable(fmm ${MY_SRC})

target_link_libraries(fmm autopas)

# add the test to checkExamples
add_dependencies(checkExamples fmm)

# -----------------test-----------------
# add check for current target
# cmake-format: off
add_test(
        NAME fmm.test
        COMMAND fmm
        CONFIGURATIONS checkExamples
)
# cmake-format: on
