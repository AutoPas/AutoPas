# needed for GCC to vectorize LJFunctor.SoAFunctor
option(ENABLE_FAST_MATH "Sets --ffast-math which is needed for gcc to vectoize efficiently" OFF)
if (ENABLE_FAST_MATH)
    message(WARNING "Fast-Math might cause particle loss! Only use this if you know what you are doing!")
endif()

target_compile_options(autopas
        PUBLIC
        # Needed to vectorize sqrt()
        -fno-math-errno
        # fast math for better vectorization
        $<$<AND:$<BOOL:${ENABLE_FAST_MATH}>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-ffast-math>
        # INTEL: per default fast math is on. Disable via fp-model precise
        $<$<AND:$<NOT:$<BOOL:${ENABLE_FAST_MATH}>>,$<CXX_COMPILER_ID:Intel>>:-fp-model precise>
        # Warnings:
        # no warnings for intel because it's mainly spam
        $<$<CXX_COMPILER_ID:GNU>:-Wsuggest-override -Wall -Wno-unused-variable -Wno-unused-function>
        $<$<CXX_COMPILER_ID:Clang>:-Wall>
        # @TODO clean up code with -Weffc++
        )

message(STATUS "fno-math-errno set. This is needed to vectorize, e.g., sqrt().")
