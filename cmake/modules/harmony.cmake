message(STATUS "harmony - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
    harmony
    URL
        ${CMAKE_SOURCE_DIR}/libs/harmony.zip
    URL_HASH MD5=b4b6ff48da4c509276a0f89e5656dab8
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/harmony/include/lib/libharmony.a
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/harmony
    # since we only unpack a header lib src == include
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/harmony/include
    # tell cmake to run make inside the source folder
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(harmony source_dir)

add_dependencies(autopas harmony)

target_include_directories(autopas SYSTEM PUBLIC ${source_dir})

target_compile_definitions(autopas PUBLIC HARMONY_HOME="HARMONY_HOME=${CMAKE_CURRENT_BINARY_DIR}/harmony/include")