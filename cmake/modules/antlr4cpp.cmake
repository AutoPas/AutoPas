message(STATUS "antlr4cpp - using bundled version")

include(ExternalProject)

# check if uuid is installed on the system, since uuid-dev is a dependency of antlr4cpp
find_package(PkgConfig)
pkg_check_modules(UUID QUIET uuid)

# if uuid was not found on system we download and install it locally
if (NOT UUID_FOUND)
    message(STATUS "UUID not found - using bundled version")

    set(LIBUUID_INSTALL_DIR "${CMAKE_BINARY_DIR}/uuid/install")
    set(UUID_CONFIG_PARAMS "--prefix=${LIBUUID_INSTALL_DIR}")
    ExternalProject_Add(
        uuid_bundled
        PREFIX uuid
        URL https://deac-ams.dl.sourceforge.net/project/libuuid/libuuid-1.0.3.tar.gz
        URL_HASH MD5=d44d866d06286c08ba0846aba1086d68
        BUILD_IN_SOURCE TRUE
        INSTALL_DIR "install"
        CONFIGURE_COMMAND "./configure" ${UUID_CONFIG_PARAMS}
        BUILD_COMMAND "make"
    )
    
else()
    set(LIBUUID_INSTALL_DIR "")
    message(STATUS "UUID found - using system version")
    # add a dummy target so the dependency in antlr4cpp_bundled is fulfilled if uuid was found on the system
    add_custom_target(uuid_bundled)
endif ()

ExternalProject_ADD(
        antlr4cpp_bundled
        PREFIX           antlr4cppPrefix
        URL              ${PROJECT_SOURCE_DIR}/libs/antlr4-cpp-runtime-4.9.3-source.zip
        URL_HASH         MD5=eafa4fef583e12e963062882773461be
        # pass PKG_CONFIG_PATH as a environment variable to cmake so find_package() in antlr4cpp's CMakeLists.txt can find uuid if using the bundled version
        CMAKE_COMMAND    PKG_CONFIG_PATH=${LIBUUID_INSTALL_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH} cmake
        BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/antlr4cppPrefix/install/lib/libantlr4-runtime.a
        CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/antlr4cppPrefix/install -DCMAKE_CXX_FLAGS=-w
        # pass LIBRARY_PATH as a environment variable to make so the linker finds the uuid-lib if using the bundled version
        BUILD_COMMAND    LIBRARY_PATH=${LIBUUID_INSTALL_DIR}/lib:$ENV{LIBRARY_PATH} make
        DEPENDS          uuid_bundled
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