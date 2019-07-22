file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
    "*.cuh"
    "*.cu"
    "../md-flexible/ParticleClassLibrary.h"
    "../md-flexible/ParticleClassLibrary.cpp"
)

if (NOT ${AUTOPAS_ENABLE_CUDA})
    message(STATUS "cuda disabled, not building bench-cuda")
    return()
endif ()

add_executable(bench-cuda ${MY_SRC})

target_link_libraries(bench-cuda autopas)

# -----------------test-----------------
# add check for current target
add_test(NAME bench-cuda.test COMMAND bench-cuda CONFIGURATIONS checkExamples)

# add the test to checkExamples
add_dependencies(checkExamples bench-cuda)
