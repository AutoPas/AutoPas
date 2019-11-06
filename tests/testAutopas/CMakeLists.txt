file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
    "../../examples/md-flexible/parsing/MDFlexConfig.cpp"
    "../../examples/md-flexible/parsing/MDFlexConfig.h"
    "../../examples/md-flexible/parsing/YamlParser.cpp"
    "../../examples/md-flexible/parsing/YamlParser.h"
    "../../examples/md-flexible/Objects/Objects.h"
    "../../examples/md-flexible/Generator.h"
    "../../examples/md-flexible/Thermostat.h"
    "../../examples/md-flexible/TimeDiscretization.h"
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

target_compile_definitions(
    runTests PRIVATE
    YAMLDIRECTORY=\"${PROJECT_SOURCE_DIR}/tests/testAutopas/exampleTests/md-flexible/yamlTestFiles/\"
)

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
