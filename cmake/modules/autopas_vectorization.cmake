option(
    AUTOPAS_USE_AUTOVEC "Enable generations of SIMD vector instructions" ON
)
if (AUTOPAS_USE_AUTOVEC)
    message(STATUS "Vectorization enabled.")

    target_compile_options(
        autopas
        PUBLIC
            # OpenMP SIMD pragmas
            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-fopenmp-simd>
            $<$<OR:$<CXX_COMPILER_ID:Intel>,$<CXX_COMPILER_ID:IntelLLVM>>:-qopenmp-simd>

            # Always use native architecture
            -march=native
    )
else ()
    message(STATUS "Vectorization disabled.")

    target_compile_options(
        autopas
        PUBLIC
            $<$<CXX_COMPILER_ID:GNU>:-fno-tree-vectorize>
            $<$<CXX_COMPILER_ID:Clang>:-fno-vectorize>
            $<$<CXX_COMPILER_ID:Intel>:-no-vec>
    )
endif ()
