FILE(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )
if(ENABLE_CUDA)
FILE(GLOB_RECURSE CU_SRC
        "*.cu"
        "*.cuh"
        )
endif()

ADD_EXECUTABLE(runTests
        ${MY_SRC}
        $<$<BOOL:${ENABLE_CUDA}>:${CU_SRC}>
        )

TARGET_LINK_LIBRARIES(runTests
        autopas
        libgtest
        libgmock
        )

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
GTEST_DISCOVER_TESTS(runTests)
