file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        "*.cuh"
        "*.cu"
        )

if (NOT ${ENABLE_CUDA})
    message(STATUS "cuda disabled, not building bench-cuda")
    return()
endif()

ADD_EXECUTABLE(bench-cuda
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(bench-cuda
        autopas)

#-----------------test-----------------
# add check for current target
add_test(
        NAME bench-cuda.test
        COMMAND bench-cuda
        CONFIGURATIONS checkExamples
)

# add the test to checkExamples
add_dependencies(checkExamples bench-cuda)
