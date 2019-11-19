# first try: check if we find any installed version
find_package(spdlog QUIET)
if (spdlog_FOUND AND "${spdlog_VERSION}" VERSION_GREATER_EQUAL 1.3.1)
    message(STATUS "spdlog - using installed system version ${spdlog_VERSION}")
    target_link_libraries(autopas PUBLIC spdlog::spdlog)
else ()
    # system version not found -> install bundled version
    message(STATUS "spdlog - not found or version older than 1.3.1")
    message(
        STATUS
            "spdlog - if you want to use your version point the cmake variable spdlog_DIR to the directory containing spdlogConfig.cmake in order to pass hints to find_package"
    )
    message(STATUS "spdlog - using bundled version 1.4.3 (commit 79259fd)")

    include(ExternalProject)
    ExternalProject_Add(
        spdlog_external
        URL
            # spdlog master:
            # https://github.com/gabime/spdlog/archive/v1.x.zip
            # spdlog commit 79259fd:
            ${CMAKE_SOURCE_DIR}/libs/spdlog-1.x.zip
        URL_HASH MD5=7415a9768f3433bd93d78c1c87fd0576
        BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/spdlog/src/spdlog-build/libspdlog.a
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/spdlog
        # Disable stuff we don't need. Especially warnings.
        CMAKE_ARGS
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_BUILD_TYPE=RELEASE
            -DSPDLOG_BUILD_EXAMPLE=OFF
            -DSPDLOG_BUILD_TESTS=OFF
            -DCMAKE_CXX_FLAGS=-w
            # Disable install step
        INSTALL_COMMAND ""
    )
    ExternalProject_Get_Property(
        spdlog_external
        source_dir
        binary_dir
        install_dir
    )

    add_library(
        spdlog::spdlog
        STATIC
        IMPORTED
        GLOBAL
    )
    add_dependencies(spdlog::spdlog spdlog_external)

    # create directory otherwise cmake will complain during generate step bc it would only be
    # created by make
    file(MAKE_DIRECTORY "${install_dir}/src/spdlog_external/include")

    # define interesting
    set_target_properties(
        spdlog::spdlog
        PROPERTIES
            "IMPORTED_LOCATION"
            "${binary_dir}/libspdlog.a"
            "INTERFACE_INCLUDE_DIRECTORIES"
            "${install_dir}/src/spdlog_external/include"
    )

endif ()
