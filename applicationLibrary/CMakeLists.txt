file(
        GLOB_RECURSE
        MY_SRC
        "molecularDynamics/molecularDynamicsLibrary/*.h"
)

add_library(molecularDynamicsLibrary STATIC ${MY_SRC})

target_link_libraries(molecularDynamicsLibrary PUBLIC autopas)

target_include_directories(molecularDynamicsLibrary PUBLIC ./molecularDynamics)