cmake_minimum_required(VERSION 3.3)

CMAKE_POLICY(SET CMP0057 NEW)


project(autopas)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.1)
    #Wsuggest-override only exists for g++ starting at version 5.1
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
endif ()

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

ENABLE_TESTING()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(tests)

include(version.cmake)
include(clang-tidy.cmake)
include(clang-format.cmake)
include(doxygen.cmake)