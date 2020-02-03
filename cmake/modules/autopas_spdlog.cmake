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

# Enable FetchContent CMake module
include(FetchContent)

# Build spdlog and make the cmake targets available
FetchContent_Declare(
    spdlog
    URL
    # spdlog master:
    # https://github.com/gabime/spdlog/archive/v1.x.zip
    # spdlog commit 79259fd:
    ${AUTOPAS_SOURCE_DIR}/libs/spdlog-1.x.zip
    URL_HASH MD5=7415a9768f3433bd93d78c1c87fd0576
)

# Disable stuff we don't need
option(SPDLOG_BUILD_EXAMPLE "" OFF)
option(SPDLOG_BUILD_TESTS "" OFF)
# hide options from ccmake
mark_as_advanced(SPDLOG_BUILD_EXAMPLE)
mark_as_advanced(SPDLOG_BUILD_TESTS)

FetchContent_MakeAvailable(spdlog)

# Disable warnings
target_compile_options(spdlog PRIVATE -w)
