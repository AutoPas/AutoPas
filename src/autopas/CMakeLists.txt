file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_LIBRARY(autopas STATIC
        ${MY_SRC}
        )

target_link_libraries(autopas
        PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${OPENMP}>:OpenMP::OpenMP_CXX>
        )

target_compile_definitions(autopas
        PUBLIC
        $<$<BOOL:${OPENMP}>:AUTOPAS_OPENMP>
        $<$<BOOL:${VERBOSE_LOGGING}>:AUTOPAS_VERBOSE_LOG>
)

TARGET_INCLUDE_DIRECTORIES(autopas PUBLIC ../.)
