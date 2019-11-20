file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)
if (AUTOPAS_ENABLE_CUDA)
    file(
        GLOB_RECURSE
        CU_SRC
        "*.cu"
        "*.cuh"
    )
endif ()

add_executable(runTests $<$<BOOL:${AUTOPAS_ENABLE_CUDA}>:${CU_SRC}> ${MY_SRC})

target_compile_definitions(runTests PRIVATE)

target_link_libraries(
    runTests
    autopas
    libgtest
    libgmock
    libyaml-cpp
)

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
gtest_discover_tests(runTests)
