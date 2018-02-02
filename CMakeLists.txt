cmake_minimum_required(VERSION 3.3)

CMAKE_POLICY(SET CMP0057 NEW)


project(autopas)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.1)
	#Wsuggest-override only exists for g++ starting at version 5.1
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
endif()

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

ENABLE_TESTING()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(tests)



if ( CMAKE_VERSION GREATER "3.5" )
    set(ENABLE_CLANG_TIDY OFF CACHE BOOL "Add clang-tidy automatically to builds")
    #set(ENABLE_CLANG_TIDY ON BOOL "Add clang-tidy automatically to builds")
    if (ENABLE_CLANG_TIDY)
        find_program (CLANG_TIDY_EXE
                NAMES "clang-tidy"
                DOC "Path to clang-tidy executable" )
        if (CLANG_TIDY_EXE)
            message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
            set(CLANG_TIDY_CHECKS "-*,modernize-*,-clang-analyzer-alpha.*")
                message(STATUS "cmake source dir: ${CMAKE_SOURCE_DIR}")
            set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE};-checks=${CLANG_TIDY_CHECKS};-header-filter='${CMAKE_SOURCE_DIR}/*'"
                    CACHE STRING "" FORCE)
        else()
            message(AUTHOR_WARNING "clang-tidy not found!")
            set(CMAKE_CXX_CLANG_TIDY "" CACHE STRING "" FORCE) # delete it
        endif()
    else()
        message(STATUS "clang tidy disabled, to enable specify -DENABLE_CLANG_TIDY=\"ON\" when calling cmake")
    endif()
else()
    message(STATUS "cmake <= 3.5")
endif()