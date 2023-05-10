file(
        GLOB_RECURSE
        MY_SRC
        "*.cpp"
        "*.h"
)

file(
        GLOB_RECURSE
        CU_SRC
        "*.cu"
        "*.cuh"
)

add_library(molecularDynamicsLibrary STATIC ${MY_SRC})

target_link_libraries(molecularDynamicsLibrary PUBLIC autopas)

target_include_directories(molecularDynamicsLibrary PUBLIC ./)
