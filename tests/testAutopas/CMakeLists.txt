file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
        "../../examples/md-flexible/YamlParser.cpp"
        "../../examples/md-flexible/YamlParser.h"
)
if (AUTOPAS_ENABLE_CUDA)
    file(
        GLOB_RECURSE
        CU_SRC
        "*.cu"
        "*.cuh"
    )
endif ()

    configure_file(${PROJECT_SOURCE_DIR}/tests/testAutopas/exampleTests/md-flexible/testParsing.yaml ${CMAKE_BINARY_DIR}/tests/testAutopas COPYONLY)

add_executable(runTests $<$<BOOL:${AUTOPAS_ENABLE_CUDA}>:${CU_SRC}> ${MY_SRC})
include_directories(${PROJECT_SOURCE_DIR}/libs/yaml-cpp/include)

target_link_libraries(
    runTests
    autopas
    libgtest
    libgmock
    ${PROJECT_SOURCE_DIR}/libs/yaml-cpp/lib/libyaml-cpp.a
)

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
gtest_discover_tests(runTests)
