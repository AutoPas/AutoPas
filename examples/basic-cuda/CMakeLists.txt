file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        "*.cuh"
        "*.cu"
        )

if (NOT ${ENABLE_CUDA})
    message(STATUS "cuda disabled, not building basic-cuda")
    return()
endif()

ADD_EXECUTABLE(basic-cuda
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(basic-cuda
        autopas)

#-----------------test-----------------
# add check for current target
add_test(
        NAME basic-cuda.test
        COMMAND basic-cuda
        CONFIGURATIONS checkExamples
)

# add the test to checkExamples
add_dependencies(checkExamples basic-cuda)
