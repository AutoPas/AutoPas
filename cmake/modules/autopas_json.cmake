# Gets nlohmann/json, which only python-based tuning needs, by (in order of priority): reusing a target a parent
# project provides, an installed version via find_package, or the bundled 3.12.0 copy in libs/json.
# Set json_ForceBundled=ON to force bundled.
#
# Both `nlohmann_json::nlohmann_json` (canonical) and `nlohmann_json` are made to resolve to the same target.

if (AUTOPAS_ENABLE_PYTHON_BASED_TUNING)
    option(json_ForceBundled "Ignore any provided/installed nlohmann/json and always use the bundled copy." OFF)
    mark_as_advanced(json_ForceBundled)

    if (NOT json_ForceBundled AND NOT AUTOPAS_FORCE_ALL_BUNDLED)
        set(expectedVersion 3.12.0)
        # Path 1: reuse a json target a parent project already defined.
        if (TARGET nlohmann_json::nlohmann_json OR TARGET nlohmann_json)
            message(STATUS "AutoPas: Reusing nlohmann/json provided by parent project")
            autopas_warn_if_parent_version_too_old(json nlohmann_json ${expectedVersion} nlohmann_json::nlohmann_json nlohmann_json)
            autopas_alias_dependency(nlohmann_json nlohmann_json::nlohmann_json)
            return()
        endif ()
        # Path 2: installed version; the version arg enforces our minimum.
        find_package(nlohmann_json ${expectedVersion} QUIET)
        if (nlohmann_json_FOUND)
            message(STATUS "json - using installed version ${nlohmann_json_VERSION}")
            autopas_promote_global(nlohmann_json::nlohmann_json)
            autopas_alias_dependency(nlohmann_json nlohmann_json::nlohmann_json)
            return()
        endif ()
        message(STATUS "json - no installed version >= ${expectedVersion} found; using bundled copy")
    else ()
        autopas_error_if_forced_bundled_collides(json json_ForceBundled nlohmann_json::nlohmann_json nlohmann_json)
    endif ()

    # Path 3 + fallback: bundled 3.12.0
    message(STATUS "json - using bundled version 3.12.0")

    # Do not add_subdirectory, as this is a header only lib. Just define a library from the single header.
    add_library(nlohmann_json INTERFACE)

    target_include_directories(nlohmann_json SYSTEM INTERFACE "${AUTOPAS_SOURCE_DIR}/libs/json/single_include")

    # Publish the canonical nlohmann_json::nlohmann_json name alongside the plain nlohmann_json target.
    autopas_alias_dependency(nlohmann_json nlohmann_json::nlohmann_json)
endif ()