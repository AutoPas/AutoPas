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
        $<$<AND:$<BOOL:${ENABLE_FAST_MATH}>,$<CXX_COMPILER_ID:Intel>>:-fp-model precise>
        # Wsuggest-override only exists for g++ starting at version 5.1
        $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,5.1>>:-Wsuggest-override -Wall -Wno-unused-variable -Wno-unused-function>
        # @TODO clean up code with -Weffc++
        # -Weffc++
)

message(STATUS "fno-math-errno set. This is needed to vectorize, e.g., sqrt().")
