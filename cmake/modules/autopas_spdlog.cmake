option(spdlog_ForceBundled "Do not look for an installed version, always use bundled." ON)

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
        # spdlog commit e86be93 (15.11.2021):
        ${AUTOPAS_SOURCE_DIR}/libs/spdlog-1.x.zip
    URL_HASH MD5=77292ebfc86717e1b5914c4d7b69140f
)

# Disable stuff we don't need
option(SPDLOG_BUILD_EXAMPLE "" OFF)
option(SPDLOG_BUILD_TESTS "" OFF)
option(SPDLOG_INSTALL "" OFF)

# hide options from ccmake
mark_as_advanced(
    SPDLOG_BUILD_BENCH
    SPDLOG_BUILD_EXAMPLE
    SPDLOG_BUILD_EXAMPLE_HO
    SPDLOG_BUILD_SHARED
    SPDLOG_BUILD_TESTS
    SPDLOG_BUILD_TESTS_HO
    SPDLOG_CLOCK_COARSE
    SPDLOG_FMT_EXTERNAL
    SPDLOG_INSTALL
    SPDLOG_NO_ATOMIC_LEVELS
    SPDLOG_NO_EXCEPTIONS
    SPDLOG_NO_THREAD_ID
    SPDLOG_NO_TLS
    SPDLOG_PREVENT_CHILD_FD
    SPDLOG_SANITIZE_ADDRESS
)

FetchContent_MakeAvailable(spdlog)

if (IS_DIRECTORY "${spdlog_SOURCE_DIR}")
    set_property(DIRECTORY ${spdlog_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
endif ()

# let ccmake and cmake-gui offer options
set(AUTOPAS_MIN_LOG_LVL
        "INFO"
        CACHE
        STRING "Choose the finest log level to be compiled."
)
set_property(CACHE AUTOPAS_MIN_LOG_LVL PROPERTY STRINGS "TRACE;DEBUG;INFO;WARN;ERROR;CRITICAL;OFF")

# Disable warnings
target_compile_options(spdlog PRIVATE -w)
# Set the finest compiled log level on the cmake target. Everything that includes this target will be affected!
target_compile_options(spdlog PUBLIC -DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${AUTOPAS_MIN_LOG_LVL})

get_target_property(propval spdlog INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(spdlog SYSTEM PUBLIC "${propval}")
