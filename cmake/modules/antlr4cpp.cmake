message(STATUS "antlr4cpp - using bundled version")

include(ExternalProject)

ExternalProject_ADD(
        antlr4cpp_bundled
        PREFIX      antlr4cppPrefix
        URL         ${PROJECT_SOURCE_DIR}/libs/antlr4-cpp-runtime-4.9.2-source.zip
        URL_HASH    MD5=187bb8d0ecb4410fe0986a735021c4f1
        BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/antlr4cppPrefix/install/lib/libantlr4-runtime.a
        CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/antlr4cppPrefix/install -DCMAKE_CXX_FLAGS=-w
)

# create dummy target that contains all information to easily link against
add_library(antlr4cpp
        STATIC
        IMPORTED
        GLOBAL
        )

add_dependencies(antlr4cpp antlr4cpp_bundled)

ExternalProject_Get_Property(antlr4cpp_bundled install_dir)
set_target_properties(
        antlr4cpp
        PROPERTIES "IMPORTED_LOCATION" "${install_dir}/install/lib/libantlr4-runtime.a"
)

# create directory otherwise cmake will complain during generate step since this is only generated during make
file(MAKE_DIRECTORY "${install_dir}/install/include/antlr4-runtime")

target_include_directories(
        antlr4cpp SYSTEM
        INTERFACE "${install_dir}/install/include/antlr4-runtime"
)
