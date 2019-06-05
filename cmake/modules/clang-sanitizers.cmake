option(AUTOPAS_ENABLE_ADDRESS_SANITIZER "Adds clang's address sanitizer to the build options" OFF)
option(AUTOPAS_ENABLE_MEMORY_SANITIZER "Adds clang's memory sanitizer to the build options" OFF)
option(AUTOPAS_ENABLE_THREAD_SANITIZER "Adds clang's thread sanitizer to the build options" OFF)

if (
    AUTOPAS_ENABLE_ADDRESS_SANITIZER
    OR AUTOPAS_ENABLE_MEMORY_SANITIZER
    OR AUTOPAS_ENABLE_THREAD_SANITIZER
)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer")
else ()
    message(STATUS "clang sanitizers disabled")
endif ()

if (AUTOPAS_ENABLE_ADDRESS_SANITIZER)
    message(STATUS "ADDRESS SANITIZER ENABLED!!!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
else ()
    message(STATUS "ADDRESS SANITIZER DISABLED")
endif ()

if (AUTOPAS_ENABLE_MEMORY_SANITIZER)
    message(STATUS "MEMORY SANITIZER ENABLED!!!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=memory")
else ()
    message(STATUS "MEMORY SANITIZER DISABLED")
endif ()

if (AUTOPAS_ENABLE_THREAD_SANITIZER)
    message(STATUS "THREAD SANITIZER ENABLED!!!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=thread")
else ()
    message(STATUS "THREAD SANITIZER DISABLED")
endif ()
