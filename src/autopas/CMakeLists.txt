file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_LIBRARY(autopas STATIC
        ${MY_SRC}
        containers/directSum/DirectSumKokkosTraversal.h)

target_link_libraries(autopas
        PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        # Do not link against openmp when using archer
        $<$<AND:$<BOOL:${OPENMP}>,$<NOT:$<BOOL:${ARCHER}>>>:OpenMP::OpenMP_CXX>
        )



# Ompstuff needs to be here because OpenMP.cmake needs to run before this file to create to OpenMP target.
# this can be resolved by upgrading to CMake 3.13 and enforcing CMP0079.
target_compile_definitions(autopas
        PUBLIC
        $<$<BOOL:${OPENMP}>:AUTOPAS_OPENMP>
        )

TARGET_INCLUDE_DIRECTORIES(autopas PUBLIC ../.)
