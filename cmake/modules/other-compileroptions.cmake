#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.1)
    #Wsuggest-override only exists for g++ starting at version 5.1
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override -Wall -Wno-unused-variable -Wno-unused-function")
endif ()

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

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")