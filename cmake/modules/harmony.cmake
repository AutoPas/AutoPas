message(STATUS "harmony - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
    harmony
    URL
        ${CMAKE_SOURCE_DIR}/libs/harmony.zip
    URL_HASH MD5=b4b6ff48da4c509276a0f89e5656dab8
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/harmony
    # since we only unpack a header lib src == include
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/harmony/include
    # Disable build & install steps
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(harmony source_dir)

add_dependencies(autopas harmony)

target_include_directories(autopas SYSTEM PUBLIC ${source_dir})
