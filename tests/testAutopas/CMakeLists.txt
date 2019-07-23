# just for faster building:
file(
    GLOB
    MY_SRC
    "*.cpp"
    "*.h"
    "exampleTests/md-flexible/*.cpp"
    "exampleTests/md-flexible/*.h"
    "../../examples/md-flexible/ParticleClassLibrary.cpp"
    "../../examples/md-flexible/ParticleClassLibrary.h"
)

# file( GLOB_RECURSE MY_SRC "*.cpp" "*.h" "../../examples/md-flexible/ParticleClassLibrary.cpp"
# "../../examples/md-flexible/ParticleClassLibrary.h" )
if (AUTOPAS_ENABLE_CUDA)
    file(
        GLOB_RECURSE
        CU_SRC
        "*.cu"
        "*.cuh"
    )
endif ()

add_executable(runTests $<$<BOOL:${AUTOPAS_ENABLE_CUDA}>:${CU_SRC}> ${MY_SRC})

target_link_libraries(
    runTests
    autopas
    libgtest
    libgmock
)

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
gtest_discover_tests(runTests)
