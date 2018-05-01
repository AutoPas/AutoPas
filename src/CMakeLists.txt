file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_LIBRARY(autopas STATIC
        ${MY_SRC}
        )

find_package(Threads REQUIRED)

target_link_libraries(autopas
        PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        )

TARGET_INCLUDE_DIRECTORIES(autopas PUBLIC .)
