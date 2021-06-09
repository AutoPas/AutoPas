# needed for GCC to vectorize LJFunctor.SoAFunctor
option(
    AUTOPAS_ENABLE_FAST_MATH "Sets --ffast-math which is needed for gcc to vectoize efficiently" OFF
)
if (AUTOPAS_ENABLE_FAST_MATH)
    message(
        WARNING
            "Fast-Math might cause particle loss! Only use this if you know what you are doing!"
    )
endif ()

# autopas requires c++17. If cmake < 3.17 is used this is set globally in the top level
# CMakeLists.txt
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
    target_compile_features(autopas PUBLIC cxx_std_17)
endif ()

target_compile_options(
    autopas
    PUBLIC
        # Needed to vectorize sqrt()
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-fno-math-errno
        # fast math for better vectorization
        $<$<AND:$<BOOL:${AUTOPAS_ENABLE_FAST_MATH}>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-ffast-math>
        # Clang: set OpenMP version to 4.5
        $<$<CXX_COMPILER_ID:Clang>:$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-fopenmp-version=45>
        # INTEL: per default fast math is on. Disable via fp-model precise
        $<$<AND:$<NOT:$<BOOL:${AUTOPAS_ENABLE_FAST_MATH}>>,$<CXX_COMPILER_ID:Intel>>:$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-fp-model
        precise>
        # Warnings:
    PRIVATE
        # no warnings for intel because it's mainly spam, but we disable one, because of a compiler
        # bug:
        # https://software.intel.com/en-us/forums/intel-c-compiler/topic/814098
        $<$<CXX_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-wd13212>
        $<$<CXX_COMPILER_ID:GNU>:
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-Wsuggest-override
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-Wall
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-Wno-unused-variable
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-Wno-unused-function
        >
        $<$<CXX_COMPILER_ID:Clang>:
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-Wall
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-Wextra
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>-Wno-unused-parameter     # triggered by functions with disabled bodies
        >
        # @TODO clean up code with -Weffc++
)

message(STATUS "fno-math-errno set. This is needed to vectorize, e.g., sqrt().")
