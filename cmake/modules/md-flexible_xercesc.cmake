message(STATUS "xcerces-c - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
    xerces
    URL
        # xerces:
        # https://www-eu.apache.org/dist//xerces/c/3/sources/xerces-c-3.2.2.zip
        # xerces-c-3.2.2:
        ${CMAKE_SOURCE_DIR}/libs/xerces-c-3.2.2.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/xerces
    # since we only unpack a header lib src =
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/xerces/
    CONFIGURE_COMMAND
    BUILD_COMMAND
    INSTALL_COMMAND
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(xerces source_dir)
add_dependencies(md-flexible xerces)

# add_library(libyaml ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp/src/yamlcpp-build/libyaml-cpp.a)
# add_library(libyaml ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp/include)


target_include_directories(md-flexible SYSTEM PUBLIC ${source_dir})
# target_link_libraries(md-flexible ${source_dir}) target_link_libraries(md-flexible libyaml)

# set(libyaml ${source_dir}/libyaml-cpp.a) target_link_libraries(md-flexible libyaml)
