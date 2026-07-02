set(AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
        # Default is OFF just for faster default compilation time.
        OFF
        CACHE
        BOOL "Enables rules-based tuning and fuzzy tuning, which compiles the bundled ANTLR cpp runtime."
        )

if (AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING)
    message(STATUS "Rules-Based Tuning Enabled")
    message(STATUS "antlr4cpp - using bundled subtree (libs/antlr4 at tag 4.13.2, Cpp runtime only)")

    include(FetchContent)

    # Suppress antlr4cpp's tests, shared lib, and demo before pulling in its CMakeLists.
    set(ANTLR_BUILD_CPP_TESTS OFF CACHE BOOL "" FORCE)
    set(ANTLR_BUILD_SHARED    OFF CACHE BOOL "" FORCE)
    set(WITH_DEMO             OFF CACHE BOOL "" FORCE)
    set(DISABLE_WARNINGS      ON  CACHE BOOL "" FORCE)

    FetchContent_Declare(
            autopas_antlr4cpp
            SOURCE_DIR ${PROJECT_SOURCE_DIR}/libs/antlr4/runtime/Cpp
    )
    FetchContent_MakeAvailable(autopas_antlr4cpp)

    # antlr4cpp's CMakeLists produces a static library target named `antlr4_static`.
    # Alias it to `antlr4cpp` to preserve the link name used in src/autopas/CMakeLists.txt.
    add_library(antlr4cpp ALIAS antlr4_static)

    # Hide antlr4cpp's cache options from ccmake/cmake-gui.
    mark_as_advanced(
            ANTLR_BUILD_CPP_TESTS
            ANTLR_BUILD_SHARED
            ANTLR_BUILD_STATIC
            DISABLE_WARNINGS
            WITH_DEMO
            WITH_LIBCXX
            WITH_STATIC_CRT
            TRACE_ATN
    )

    # Keep antlr4cpp targets out of the default "all" build/install step.
    if (IS_DIRECTORY "${autopas_antlr4cpp_SOURCE_DIR}")
        set_property(DIRECTORY ${autopas_antlr4cpp_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
    endif ()

else ()
    message(STATUS "Rules-Based Tuning Disabled. Bundled ANTLR will not be compiled.")
endif ()