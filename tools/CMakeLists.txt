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

add_library(autopas-tools STATIC $<$<BOOL:${AUTOPAS_ENABLE_CUDA}>:${CU_SRC}> ${MY_SRC})

target_link_libraries(
        autopas-tools
        PUBLIC
        autopas
)

target_include_directories(autopas-tools PUBLIC ./)
