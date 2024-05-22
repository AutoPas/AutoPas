set(AUTOPAS_ENABLE_RULES_BASED_TUNING
        ON
        CACHE
        BOOL "Enables rules-based tuning, which, if using the bundled version, will compile ANTLR and, if uuid is not installed already, will compile uuid which can cause errors."
        )

if (AUTOPAS_ENABLE_RULES_BASED_TUNING)
    message(STATUS "Rules-Based Tuning Enabled")
    message(STATUS "antlr4cpp - using bundled version")

    include(ExternalProject)

    # check if uuid-dev is installed on the system, since this is a dependency of antlr4cpp
    find_package(PkgConfig)
    pkg_check_modules(UUID QUIET uuid)

    # if uuid-dev was not found on system we install it locally
    if (NOT UUID_FOUND)
        message(STATUS "UUID not found - using bundled version")

        set(LIBUUID_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/uuid/install")
        set(LIBUUID_PKGCONFIG_DIR ${LIBUUID_INSTALL_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})
        set(LIBUUID_LIBRARY_DIR ${LIBUUID_INSTALL_DIR}/lib:$ENV{LIBRARY_PATH})
        set(UUID_CONFIG_PARAMS "--prefix=${LIBUUID_INSTALL_DIR}")
        ExternalProject_Add(
            uuid_bundled
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/uuid
            URL ${PROJECT_SOURCE_DIR}/libs/libuuid-1.0.3.zip
            URL_HASH MD5=9fd1e87682d24d6ca22941e0af339c8a
            BUILD_IN_SOURCE TRUE
            INSTALL_DIR "install"
            CONFIGURE_COMMAND "./configure" ${UUID_CONFIG_PARAMS}
            BUILD_COMMAND ${MAKE_EXE}
        )

    else()
        message(STATUS "UUID found - using system version")
        # add a dummy target so the dependency in antlr4cpp_bundled is fulfilled if uuid-dev was found on the system
        add_custom_target(uuid_bundled)
    endif ()

    find_package(utf8cpp QUIET)

    # if utf8cpp was not found on system we install it locally
    if (NOT utf8cpp_FOUND)
        message(STATUS "utf8cpp not found - using bundled version")

        set(UTFCPP_DIR "${CMAKE_CURRENT_BINARY_DIR}/utf8cpp")
        ExternalProject_Add(
            utf8cpp_bundled
            PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/utf8cpp
            URL             ${PROJECT_SOURCE_DIR}/libs/utfcpp-3.1.1.zip
            URL_HASH        MD5=a2cf6db2ee03ccdcf5308793400acfe1
            BUILD_IN_SOURCE TRUE
            INSTALL_DIR     "install"
            CMAKE_ARGS      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${UTFCPP_DIR}/install -DUTF8_TESTS=off -DUTF8_SAMPLES=off
        )
    else()
        message(STATUS "utf8cpp found - using system version")
        # add a dummy target so the dependency in antlr4cpp_bundled is fulfilled if utf8cpp was found on the system
        add_custom_target(utf8cpp_bundled)
    endif ()

    # install prefix for antlr
    set(antlr4cpp_prefix ${CMAKE_CURRENT_BINARY_DIR}/_deps/antlr4cppPrefix)
    # location where antlr will install its static library
    set(staticLibInstallLocation ${antlr4cpp_prefix}/install/lib/libantlr4-runtime.a)

    ExternalProject_ADD(
            antlr4cpp_bundled
            PREFIX           ${antlr4cpp_prefix}
            URL              ${PROJECT_SOURCE_DIR}/libs/antlr4-cpp-runtime-4.9.3-source.zip
            URL_HASH         MD5=eafa4fef583e12e963062882773461be
            BUILD_BYPRODUCTS ${staticLibInstallLocation}
            # pass PKG_CONFIG_PATH as a environment variable to cmake so find_package() in antlr4cpp's CMakeLists.txt can find uuid-dev if using the bundled version
            CMAKE_COMMAND    ${CMAKE_COMMAND} -E env PKG_CONFIG_PATH=${LIBUUID_PKGCONFIG_DIR} ${CMAKE_COMMAND}
            BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/antlr4cpp/install/lib/libantlr4-runtime.a
            # point antlr4cpp to utf8cpp install dir
            CMAKE_ARGS       -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${antlr4cpp_prefix}/install -DCMAKE_CXX_FLAGS=-w -DCMAKE_PREFIX_PATH=${UTFCPP_DIR}/install
            # make sure UUID and UTF8CPP is installed before antlr is installed
            DEPENDS          uuid_bundled utf8cpp_bundled
            # Patch away the shared library target
            PATCH_COMMAND    sed --in-place -e /install.*shared/,+2d               runtime/CMakeLists.txt  &&
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

else()
    message(STATUS "Rules-Based Tuning Disabled. Bundled versions of ANTLR and uuid will not be compiled.")
endif()