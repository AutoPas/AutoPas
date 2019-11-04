message(STATUS "yaml-cpp - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract yaml-cpp
ExternalProject_Add(
    yaml-cpp
    URL
        # eigen-master:
        # https://github.com/jbeder/yaml-cpp/archive/master.zip
        # commit a8ba6a8:
        ${CMAKE_SOURCE_DIR}/libs/yaml-cpp-master.zip
    URL_HASH MD5=6495f4c5dcf45414e25938e53a10ce9d
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp
    # Disable install steps
    INSTALL_COMMAND
        ""
        # Disable everything we don't need and set build type to release
    CMAKE_ARGS
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=RELEASE
        -DBUILD_TESTING=OFF
        -DYAML_CPP_BUILD_TESTS=OFF
        -DYAML_CPP_BUILD_CONTRIB=OFF
        -DYAML_CPP_BUILD_TOOLS=OFF
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(yaml-cpp install_dir binary_dir)

add_library(
    libyaml-cpp
    STATIC
    IMPORTED
    GLOBAL
)

add_dependencies(libyaml-cpp yaml-cpp)

file(MAKE_DIRECTORY "${install_dir}/src/yaml-cpp/include")

# Set libgtest properties
set_target_properties(
    libyaml-cpp
    PROPERTIES
        "IMPORTED_LOCATION"
        "${binary_dir}/libyaml-cpp.a"
        "INTERFACE_INCLUDE_DIRECTORIES"
        "${install_dir}/src/yaml-cpp/include"
)
