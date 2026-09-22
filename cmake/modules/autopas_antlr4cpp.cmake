set(AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
        # Default is OFF just for faster default compilation time.
        OFF
        CACHE
        BOOL "Enables rules-based tuning and fuzzy tuning, which, if using the bundled version, will compile ANTLR."
        )

if (AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING)
    message(STATUS "Rules-Based Tuning Enabled")
    message(STATUS "antlr4cpp - using bundled version")

    # Set ANTLR build options to disable building tests, shared libraries, and demos, and to disable warnings.
    set(ANTLR_BUILD_CPP_TESTS OFF CACHE INTERNAL "")
    set(ANTLR_BUILD_SHARED OFF CACHE INTERNAL "")
    set(WITH_DEMO             OFF CACHE INTERNAL "")
    set(DISABLE_WARNINGS      ON  CACHE INTERNAL "")

    set(ANTLR_BUILD_STATIC ON CACHE INTERNAL "")

    add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/antlr4/runtime/Cpp ${CMAKE_CURRENT_BINARY_DIR}/antlr4 EXCLUDE_FROM_ALL)

    # antlr4cpp's CMakeLists produces a static library target named `antlr4_static`.
    # Alias it to `antlr4cpp` to preserve the link name used in src/autopas/CMakeLists.txt.
    add_library(antlr4cpp ALIAS antlr4_static)

    # Hide antlr4cpp's further cmake options from ccmake/cmake-gui.
    mark_as_advanced(
            WITH_LIBCXX
            WITH_STATIC_CRT
            TRACE_ATN
    )

else()
    message(STATUS "Rules-Based Tuning Disabled. Bundled version of ANTLR will not be compiled.")
endif()