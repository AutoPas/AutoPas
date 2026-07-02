option(spdlog_ForceBundled "Do not look for an installed version, always use bundled." ON)

if (NOT ${spdlog_ForceBundled})
    # first try: check if we find any installed version
    set(expectedVersion 1.17.0)
    find_package(spdlog ${expectedVersion} QUIET)
    if (spdlog_FOUND)
        message(STATUS "spdlog - using installed system version ${spdlog_VERSION}")
        # promote target to global visibility to be used elsewhere
        set_target_properties(spdlog::spdlog PROPERTIES IMPORTED_GLOBAL TRUE)
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
message(STATUS "spdlog - using bundled version 1.17.0 (commit 79524dd)")

# Enable FetchContent CMake module
include(FetchContent)

# Build spdlog and make the cmake targets available
FetchContent_Declare(
    spdlog
    URL
        # spdlog master:
        # https://github.com/gabime/spdlog/archive/refs/tags/v1.17.0.zip
        # spdlog commit 79524dd (04.01.2026):
        ${AUTOPAS_SOURCE_DIR}/libs/spdlog-1.17.0.zip
    URL_HASH MD5=d38d278383b768847ccc4616879df42f
)

# Disable stuff we don't need (Sets values to OFF and hides them)
set(SPDLOG_BUILD_EXAMPLE OFF CACHE INTERNAL "Disable spdlog examples")
set(SPDLOG_BUILD_TESTS   OFF CACHE INTERNAL "Disable spdlog tests")
set(SPDLOG_BUILD_BENCH   OFF CACHE INTERNAL "Disable spdlog benchmarks")
set(SPDLOG_INSTALL       OFF CACHE INTERNAL "Disable spdlog install rules")
set(SPDLOG_FMT_EXTERNAL  OFF CACHE INTERNAL "Ensure bundled fmt is used")

FetchContent_MakeAvailable(spdlog)

# Hide all remaining internal SPDLOG cache variables from ccmake/cmake-gui
get_cmake_property(_all_cache_vars CACHE_VARIABLES)
foreach(_var IN LISTS _all_cache_vars)
    if(_var MATCHES "^SPDLOG_")
        mark_as_advanced(${_var})
    endif()
endforeach()

# Prevent spdlog targets from being included in the default "all" build/install step
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

# Treat spdlog headers as system headers to suppress compiler warnings originating from them
get_target_property(propval spdlog INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(spdlog SYSTEM PUBLIC "${propval}")
