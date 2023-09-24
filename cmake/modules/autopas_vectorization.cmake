option(
    AUTOPAS_USE_VECTORIZATION "Enable generations of SIMD vector instructions through omp-simd" ON
)
if (AUTOPAS_USE_VECTORIZATION)
    message(STATUS "Vectorization enabled.")
    # list of available options
    set(VECTOR_INSTRUCTIONS_OPTIONS "NATIVE;DEFAULT;SSE;AVX;AVX2;KNL")
    # set instruction set type
    set(
        AUTOPAS_VECTOR_INSTRUCTIONS
        "NATIVE"
        CACHE STRING "Vector instruction set to use (${VECTOR_INSTRUCTIONS_OPTIONS})."
    )
    # let ccmake and cmake-gui offer the options
    set_property(CACHE AUTOPAS_VECTOR_INSTRUCTIONS PROPERTY STRINGS ${VECTOR_INSTRUCTIONS_OPTIONS})

    target_compile_options(
        autopas
        PUBLIC
            # openmp simd
            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-fopenmp-simd>
            $<$<OR:$<CXX_COMPILER_ID:Intel>,$<CXX_COMPILER_ID:IntelLLVM>>:-qopenmp-simd>
            # vector instruction set
            $<$<STREQUAL:${AUTOPAS_VECTOR_INSTRUCTIONS},NATIVE>:-march=native>
            $<$<STREQUAL:${AUTOPAS_VECTOR_INSTRUCTIONS},SSE>:-msse3>
            $<$<STREQUAL:${AUTOPAS_VECTOR_INSTRUCTIONS},AVX>:-mavx>
            $<$<AND:$<STREQUAL:${AUTOPAS_VECTOR_INSTRUCTIONS},AVX2>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-mavx2
            -mfma>
            $<$<AND:$<STREQUAL:${AUTOPAS_VECTOR_INSTRUCTIONS},AVX2>,$<OR:$<CXX_COMPILER_ID:Intel>,$<CXX_COMPILER_ID:IntelLLVM>>>:-march=core-avx2
            -fma>
            $<$<AND:$<STREQUAL:${AUTOPAS_VECTOR_INSTRUCTIONS},KNL>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:IntelLLVM>>>:-march=knl>
            $<$<AND:$<STREQUAL:${AUTOPAS_VECTOR_INSTRUCTIONS},KNL>,$<CXX_COMPILER_ID:Intel>>:-xMIC-AVX512>
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
