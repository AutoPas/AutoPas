file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

add_executable(sph-main ${MY_SRC})

target_link_libraries(sph-main autopas SPHLibrary)

# -----------------test-----------------
# add check for current target
# cmake-format: off
add_test(
    NAME sph-main.test
    COMMAND sph-main | grep "time step" | tail -2 | head -1 | grep -q "time step 8"
    CONFIGURATIONS checkExamples
)
# cmake-format: on

# add the executable to checkExamples as dependency
add_dependencies(checkExamples sph-main)
