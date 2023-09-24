option(FORCE_COLORED_OUTPUT "Always produce ANSI-colored compiler output (GNU/Clang only)." FALSE)
if (${FORCE_COLORED_OUTPUT})
    target_compile_options(
        autopas
        PUBLIC
            $<$<CXX_COMPILER_ID:GNU>:-fdiagnostics-color=always>
            $<$<CXX_COMPILER_ID:Clang>:-fcolor-diagnostics>
    )
endif ()
