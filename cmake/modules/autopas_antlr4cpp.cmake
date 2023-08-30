message(STATUS "antlr4cpp - using bundled version")

include(ExternalProject)

# install prefix for antlr
set(prefix ${CMAKE_CURRENT_BINARY_DIR}/_deps/antlr4cppPrefix)
# location where antlr will install its static library
set(staticLibInstallLocation ${prefix}/install/lib/libantlr4-runtime.a)

ExternalProject_ADD(
        antlr4cpp_bundled
        PREFIX      ${prefix}
        URL         ${PROJECT_SOURCE_DIR}/libs/antlr4-cpp-runtime-4.9.3-source.zip
        URL_HASH    MD5=eafa4fef583e12e963062882773461be
        BUILD_BYPRODUCTS ${staticLibInstallLocation}
        CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${prefix}/install -DCMAKE_CXX_FLAGS=-w
        # Patch away the shared library target
        PATCH_COMMAND sed --in-place -e /install.*shared/,+2d               runtime/CMakeLists.txt  &&
                      sed --in-place -e /set_target_properties.*shared/,+9d runtime/CMakeLists.txt  &&
                      sed --in-place -e /shared/d                           runtime/CMakeLists.txt
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
        PROPERTIES "IMPORTED_LOCATION" "${staticLibInstallLocation}"
)

# create directory otherwise cmake will complain during generate step since this is only generated during make
file(MAKE_DIRECTORY "${install_dir}/install/include/antlr4-runtime")

target_include_directories(
        antlr4cpp SYSTEM
        INTERFACE "${install_dir}/install/include/antlr4-runtime"
)
