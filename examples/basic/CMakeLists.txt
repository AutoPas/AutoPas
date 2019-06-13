file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

add_executable(basic ${MY_SRC})

target_link_libraries(basic autopas)

# -----------------test-----------------
# add check for current target
add_test(NAME basic.test COMMAND basic CONFIGURATIONS checkExamples)

# add the test to checkExamples
add_dependencies(checkExamples basic)
