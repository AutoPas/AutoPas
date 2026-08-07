set(AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
        # Default is OFF just for faster default compilation time.
        OFF
        CACHE
        BOOL "Enables rules-based tuning and fuzzy tuning, which, if using the bundled version, will compile ANTLR."
        )

if (NOT AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING)
    return()
endif ()

message(STATUS "Rules-Based Tuning Enabled")

# Gets the antlr4 C++ runtime by (in order of priority): reusing a target a parent project provides, an
# installed version via find_package, or the bundled 4.13.2 copy in libs/antlr4.
# Set antlr4cpp_ForceBundled=ON to force bundled.
option(antlr4cpp_ForceBundled "Ignore any provided/installed antlr4cpp and always use the bundled copy." OFF)
mark_as_advanced(antlr4cpp_ForceBundled)

if (NOT antlr4cpp_ForceBundled AND NOT AUTOPAS_FORCE_ALL_BUNDLED)
    # Path 1: reuse an antlr4 runtime a parent project already defined, under either AutoPas's link name
    # or one of the names upstream uses.
    if (TARGET antlr4cpp OR TARGET antlr4_shared OR TARGET antlr4_static)
        message(STATUS "AutoPas: Reusing antlr4cpp provided by parent project")
        set(antlr4cppResolved TRUE)
    else ()
        # Path 2: installed version; the runtime installs a CMake package `antlr4-runtime` which defines
        # the imported targets. The version arg enforces our minimum.
        set(expectedVersion 4.13.2)
        find_package(antlr4-runtime ${expectedVersion} QUIET)
        if (antlr4-runtime_FOUND)
            message(STATUS "antlr4cpp - using installed version ${ANTLR_VERSION}")
            autopas_promote_global(antlr4_shared)
            autopas_promote_global(antlr4_static)
            set(antlr4cppResolved TRUE)
        else ()
            message(STATUS "antlr4cpp - no installed version >= ${expectedVersion} found; using bundled copy")
        endif ()
    endif ()
endif ()

if (NOT antlr4cppResolved)
    # Path 3 + fallback: bundled version.
    message(STATUS "antlr4cpp - using bundled version 4.13.2")

    # Set ANTLR build options to disable building tests, shared libraries, and demos, and to disable warnings.
    set(ANTLR_BUILD_CPP_TESTS OFF CACHE INTERNAL "")
    set(ANTLR_BUILD_SHARED OFF CACHE INTERNAL "")
    set(WITH_DEMO             OFF CACHE INTERNAL "")
    set(DISABLE_WARNINGS      ON  CACHE INTERNAL "")

    set(ANTLR_BUILD_STATIC ON CACHE INTERNAL "")

    add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/antlr4/runtime/Cpp ${CMAKE_CURRENT_BINARY_DIR}/antlr4 EXCLUDE_FROM_ALL)

    # Hide antlr4cpp's further cmake options from ccmake/cmake-gui.
    mark_as_advanced(
            WITH_LIBCXX
            WITH_STATIC_CRT
            TRACE_ATN
    )
endif ()

# AutoPas links the runtime as `antlr4cpp` (see src/autopas/CMakeLists.txt), so expose it under that name, but 
# antlr4cpp's CMakeLists produces a shared and/or static library, `antlr4_shared` and `antlr4_static`. By default, we
# alias the static library as `antlr4cpp`, but a parent project or installed package may only provide the shared
# library, in which case we alias that as `antlr4cpp`.
if (NOT TARGET antlr4cpp)
    if (TARGET antlr4_static)
        autopas_alias_dependency(antlr4cpp antlr4_static)
    else ()
        autopas_alias_dependency(antlr4cpp antlr4_shared)
    endif ()
endif ()