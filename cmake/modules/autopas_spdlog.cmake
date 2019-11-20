option(spdlog_ForceBundled "Do not look for an installed version, always use bundled." OFF)

if (NOT ${spdlog_ForceBundled})
    # first try: check if we find any installed version
    set(expectedVersion 1.3.1)
    find_package(spdlog ${expectedVersion} QUIET)
    if (spdlog_FOUND)
        message(STATUS "spdlog - using installed system version ${spdlog_VERSION}")
        # promote target to global visibility to be used elsewhere
        set_target_properties(spdlog::spdlog PROPERTIES "IMPORTED_GLOBAL" "TRUE")
        return()
    else ()
            message(STATUS "spdlog - no system version compatible to version ${expectedVersion} found")
            message(
                STATUS
                    "spdlog - if you want to use your version point the cmake variable spdlog_DIR to the directory containing spdlogConfig.cmake in order to pass hints to find_package"
            )
    endif ()
endif ()

# system version not found -> install bundled version
message(STATUS "spdlog - using bundled version 1.4.3 (commit 79259fd)")

include(ExternalProject)
ExternalProject_Add(
    spdlog_bundled
    URL
        # spdlog master:
        # https://github.com/gabime/spdlog/archive/v1.x.zip
        # spdlog commit 79259fd:
        ${CMAKE_SOURCE_DIR}/libs/spdlog-1.x.zip
    URL_HASH MD5=7415a9768f3433bd93d78c1c87fd0576
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/spdlog/src/spdlog_bundled-build/libspdlog.a
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
    spdlog_bundled
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
add_dependencies(spdlog::spdlog spdlog_bundled)

# create directory otherwise cmake will complain during generate step bc it would only be created by
# make
file(MAKE_DIRECTORY "${install_dir}/src/spdlog_bundled/include")

# define interesting
set_target_properties(
    spdlog::spdlog
    PROPERTIES
        "IMPORTED_LOCATION"
        "${binary_dir}/libspdlog.a"
        "INTERFACE_INCLUDE_DIRECTORIES"
        "${install_dir}/src/spdlog_bundled/include"
)
