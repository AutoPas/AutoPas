file(GLOB_RECURSE
    JULIA_SRC
    "src/*.cpp"
    "src/*.h"
)
add_library(julia_bindings SHARED ${JULIA_SRC})

set(julia_bindings "/mnt/c/Users/laura/Documents/test_ordner/libcxxwrap-julia/build/")

list(APPEND CMAKE_PREFIX_PATH "/mnt/c/Users/laura/Documents/test_ordner/libcxxwrap-julia/build/")

find_package(JlCxx REQUIRED)

target_include_directories(julia_bindings PUBLIC "/mnt/c/Users/laura/Documents/test_ordner/libcxxwrap-julia/build/" JlCxx::JlCxx)

target_link_libraries(julia_bindings PUBLIC autopas JlCxx::cxxwrap_julia)
target_compile_options(julia_bindings PUBLIC -fPIC)