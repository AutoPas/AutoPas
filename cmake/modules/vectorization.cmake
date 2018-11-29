option(USE_VECTORIZATION "Enable generations of SIMD vector instructions through omp-simd" ON)
if (USE_VECTORIZATION)
    MESSAGE(STATUS "Vectorization enabled.")
    # list of available options
    set(VECTOR_INSTRUCTIONS_OPTIONS "NATIVE;SSE;AVX;AVX2;KNL")
    # set instruction set type
    set(VECTOR_INSTRUCTIONS "NATIVE" CACHE STRING "Vector instruction set to use (${VECTOR_INSTRUCTIONS_OPTIONS}).")
    # let ccmake and cmake-gui offer the options
    set_property(CACHE VECTOR_INSTRUCTIONS PROPERTY STRINGS ${VECTOR_INSTRUCTIONS_OPTIONS})

    target_compile_options(autopas
            PUBLIC
            # openmp simd
            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-fopenmp-simd>
            $<$<CXX_COMPILER_ID:Intel>:-qopenmp-simd>

            # vector instruction set
            $<$<STREQUAL:${VECTOR_INSTRUCTIONS},NATIVE>:-march=native>
            $<$<STREQUAL:${VECTOR_INSTRUCTIONS},SSE>:-msse3>
            $<$<STREQUAL:${VECTOR_INSTRUCTIONS},AVX>:-mavx>
            $<$<AND:$<STREQUAL:${VECTOR_INSTRUCTIONS},AVX2>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-mavx2 -mfma>
            $<$<AND:$<STREQUAL:${VECTOR_INSTRUCTIONS},AVX2>,$<CXX_COMPILER_ID:Intel>>:-march=core-avx2 -fma>
            $<$<AND:$<STREQUAL:${VECTOR_INSTRUCTIONS},KNL>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-march=knl>
            $<$<AND:$<STREQUAL:${VECTOR_INSTRUCTIONS},KNL>,$<CXX_COMPILER_ID:Intel>>:-xMIC-AVX512>
            )

else()
	MESSAGE(STATUS "Vectorization disabled.")

    target_compile_options(autopas
            PUBLIC
            $<$<CXX_COMPILER_ID:GNU>:-fno-tree-vectorize>
            $<$<CXX_COMPILER_ID:Clang>:-fno-vectorize>
            $<$<CXX_COMPILER_ID:Intel>:-no-vec>
            )
endif ()