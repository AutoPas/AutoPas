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

add_library(autopasTools STATIC ${MY_SRC})

target_link_libraries(autopasTools PUBLIC autopas)

target_include_directories(autopasTools PUBLIC ./)
