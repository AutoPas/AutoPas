file(
        GLOB_RECURSE
        MY_SRC
        "*.cpp"
        "*.h"
)

add_executable(kokkos-test ${MY_SRC})

target_link_libraries(kokkos-test autopas)
# -----------------test-----------------
# add check for current target
#add_test(NAME basic.test COMMAND basic CONFIGURATIONS checkExamples)

# add the test to checkExamples
#add_dependencies(checkExamples basic)
