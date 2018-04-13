cmake_minimum_required(VERSION 3.3)

CMAKE_POLICY(SET CMP0057 NEW)

project(AutoPas)

set(CMAKE_CXX_STANDARD 11)
# cmake < 3.6 ignores "set(CMAKE_CXX_STANDARD 11)" when using icpc
if (CMAKE_VERSION VERSION_LESS "3.6" AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.1)
    #Wsuggest-override only exists for g++ starting at version 5.1
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
endif ()

# TODO: move this to a nicer place (own file?)
# needed for GCC to vectorize LJFunctor.SoAFunctor
option(ENABLE_FAST_MATH "Sets --ffast-math which is needed for gcc to vectoize efficiently" OFF)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (ENABLE_FAST_MATH)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if (NOT ENABLE_FAST_MATH)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model precise")
    endif ()
endif ()

################################### START VECTORIZATION ###################################
option(USE_VECTORIZATION "Enable generations of SIMD vector instructions through omp-simd" ON)
if (USE_VECTORIZATION)
    # list of available options
    set(VECTOR_INSTRUCTIONS_OPTIONS "NATIVE;SSE;AVX;AVX2;KNL")
    # set instruction set type
    set(VECTOR_INSTRUCTIONS "NATIVE" CACHE STRING "Vector instruction set to use\
 (${VECTOR_INSTRUCTIONS_OPTIONS}).")
    # let ccmake and cmake-gui offer the options
    set_property(CACHE VECTOR_INSTRUCTIONS PROPERTY STRINGS ${VECTOR_INSTRUCTIONS_OPTIONS})

    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd")

        if (VECTOR_INSTRUCTIONS MATCHES "^NATIVE$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^SSE$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^AVX$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^AVX2$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^KNL$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=knl")
        else ()
            message(SEND_ERROR "\"${VECTOR_INSTRUCTIONS}\" is an unknown vector instruction set option.\
     Available options: ${VECTOR_INSTRUCTIONS_OPTIONS}")
        endif ()

    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp-simd")

        if (VECTOR_INSTRUCTIONS MATCHES "^NATIVE$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^SSE$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^AVX$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^AVX2$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=core-avx2 -fma")
        elseif (VECTOR_INSTRUCTIONS MATCHES "^KNL$")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xMIC-AVX512")
        else ()
            message(SEND_ERROR "\"${VECTOR_INSTRUCTIONS}\" is an unknown vector instruction set option.\
     Available options: ${VECTOR_INSTRUCTIONS_OPTIONS}")
        endif ()
    endif ()
elseif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-tree-vectorize")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-vec")
    endif ()
endif ()
###################################  END VECTORIZATION  ###################################

ENABLE_TESTING()

include(address-sanitizer.cmake)
include(version.cmake)
include(clang-tidy.cmake)
include(clang-format.cmake)
include(doxygen.cmake)


ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(tests)