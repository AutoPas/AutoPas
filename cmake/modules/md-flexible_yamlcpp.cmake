message(STATUS "yamlcpp - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
        yamlcpp
        URL
        # yamlcpp:
        #GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
        #GIT_TAG yaml-cpp-0.6.2
        ${CMAKE_SOURCE_DIR}/libs/yaml-cpp.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp
        # since we only unpack a header lib src =
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp/include
        CONFIGURE_COMMAND
        BUILD_COMMAND
        INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(yamlcpp source_dir)
add_dependencies(md-flexible yamlcpp)

#add_library(libyaml ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp/src/yamlcpp-build/libyaml-cpp.a)
#add_library(libyaml ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp/include)

if(NOT "${CMAKE_CURRENT_BINARY_DIR}/yamlcpp/src/yamlcpp-build/libyaml-cpp.a")
    message(STATUS "libyaml-cpp.a FOUND")
        else()
    message(STAUS "libyaml-cpp.a NOTFOUND")
        endif()


target_include_directories(md-flexible SYSTEM PUBLIC ${source_dir})
#target_link_libraries(md-flexible ${source_dir})
#target_link_libraries(md-flexible libyaml)



#set(libyaml ${source_dir}/libyaml-cpp.a)
#target_link_libraries(md-flexible libyaml)
