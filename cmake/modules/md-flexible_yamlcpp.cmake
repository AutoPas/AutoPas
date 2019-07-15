message(STATUS "yamlcpp - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
        yamlcpp
        URL
        # yamlcpp:
        ${CMAKE_SOURCE_DIR}/libs/yaml-cpp.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp
        # since we only unpack a header lib src == include
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp/include
        # Disable build & install steps
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(yamlcpp source_dir)

add_dependencies(md-flexible yamlcpp)

target_include_directories(autopas SYSTEM PUBLIC ${source_dir})
