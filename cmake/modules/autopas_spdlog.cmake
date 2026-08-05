# Gets spdlog by (in order of priority): reusing a target a parent project provides, an installed
# version via find_package, or the bundled 1.17.0 copy in libs/spdlog. Set spdlog_ForceBundled=ON to force bundled.
option(spdlog_ForceBundled "Ignore any provided/installed spdlog and always use the bundled copy." OFF)
mark_as_advanced(spdlog_ForceBundled)

set(AUTOPAS_MIN_LOG_LVL
        "INFO"
        CACHE
        STRING "Choose the finest log level to be compiled."
)
set_property(CACHE AUTOPAS_MIN_LOG_LVL PROPERTY STRINGS "TRACE;DEBUG;INFO;WARN;ERROR;CRITICAL;OFF")

if (NOT spdlog_ForceBundled)
    # Path 1: reuse a spdlog target a parent project already defined. AUTOPAS_MIN_LOG_LVL is not applied
    # here, as the target belongs to the parent project, which decides how its spdlog is compiled.
    if (TARGET spdlog::spdlog OR TARGET spdlog)
        message(STATUS "AutoPas: Reusing spdlog provided by parent project")
        autopas_alias_dependency(spdlog spdlog::spdlog)
        return()
    endif ()
    # Path 2: installed version
    set(expectedVersion 1.17.0)
    find_package(spdlog ${expectedVersion} QUIET)
    if (spdlog_FOUND)
        message(STATUS "spdlog - using installed system version ${spdlog_VERSION}")
        autopas_promote_global(spdlog)
        autopas_promote_global(spdlog::spdlog)
        autopas_alias_dependency(spdlog spdlog::spdlog)
        # System targets can only take INTERFACE options, and are already treated as SYSTEM headers by CMake
        target_compile_options(spdlog INTERFACE -DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${AUTOPAS_MIN_LOG_LVL})
        return()
    endif ()
    message(STATUS "spdlog - no system version >= ${expectedVersion} found; using bundled copy")
endif ()

# Path 3 + fallback: bundled version.
message(STATUS "spdlog - using bundled version 1.17.0 (commit 79524dd)")

# Disable stuff we don't need (Sets values to OFF and hides them)
set(SPDLOG_BUILD_EXAMPLE OFF CACHE INTERNAL "Disable spdlog examples")
set(SPDLOG_BUILD_TESTS   OFF CACHE INTERNAL "Disable spdlog tests")
set(SPDLOG_BUILD_BENCH   OFF CACHE INTERNAL "Disable spdlog benchmarks")
set(SPDLOG_INSTALL       OFF CACHE INTERNAL "Disable spdlog install rules")
set(SPDLOG_FMT_EXTERNAL  OFF CACHE INTERNAL "Ensure bundled fmt is used")

add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/spdlog ${CMAKE_CURRENT_BINARY_DIR}/spdlog EXCLUDE_FROM_ALL)

# Hide all remaining internal SPDLOG cache variables from ccmake/cmake-gui
get_cmake_property(_all_cache_vars CACHE_VARIABLES)
foreach(_var IN LISTS _all_cache_vars)
    if(_var MATCHES "^SPDLOG_")
        mark_as_advanced(${_var})
    endif()
endforeach()

# Disable warnings
target_compile_options(spdlog PRIVATE -w)
# Set the finest compiled log level on the cmake target. Everything that includes this target will be affected!
target_compile_options(spdlog PUBLIC -DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${AUTOPAS_MIN_LOG_LVL})

# Treat spdlog headers as system headers to suppress compiler warnings originating from them
get_target_property(propval spdlog INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(spdlog SYSTEM PUBLIC "${propval}")