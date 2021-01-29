if(NOT MSVC)
    return()
endif()
PROJECT(win_getopt C)

message(STATUS "adding win_getopt library")
add_library(win_getopt OBJECT getopt.c getopt.h)
target_include_directories(win_getopt PUBLIC .)

SET_TARGET_PROPERTIES(win_getopt PROPERTIES LINKER_LANGUAGE C)