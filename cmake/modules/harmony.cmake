message(STATUS "harmony - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
    harmony_bundled
    URL
        ${CMAKE_SOURCE_DIR}/libs/harmony.zip
    URL_HASH MD5=a8768c2886bdc2e44e3b6b7d4f94729c
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/harmony/include/lib/libharmony.a
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/harmony
    # since we only unpack a header lib src == include
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/harmony/include
    # tell cmake to run make inside the source folder
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(harmony_bundled source_dir)

add_library(
        harmony
        STATIC
        IMPORTED GLOBAL
)

#target_link_libraries(harmony INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/harmony/include/lib/libharmony.a)

add_dependencies(harmony harmony_bundled)

# create directory otherwise cmake will complain during generate step bc it would only be created by make
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/harmony/include/include")

set_target_properties(
        harmony
        PROPERTIES
        "IMPORTED_LOCATION"
        "${CMAKE_CURRENT_BINARY_DIR}/harmony/include/lib/libharmony.a"
        "INTERFACE_INCLUDE_DIRECTORIES"
        "${CMAKE_CURRENT_BINARY_DIR}/harmony/include/include"
)

# Set macro needed to set environment variable for ActiveHarmony
target_compile_definitions(harmony INTERFACE HARMONY_HOME="HARMONY_HOME=${CMAKE_CURRENT_BINARY_DIR}/harmony/include")
