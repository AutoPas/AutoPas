option(OPENMP "Activates OpenMP shared memory parallelization." OFF)

if (CMAKE_CXX_COMPILER MATCHES "archer")
    message(STATUS "archer detected, OpenMP enabled by default, so skipping OpenMP package search")
    set(OPENMP ON)
    return()
endif()

if (OPENMP)
    message(STATUS "OpenMP enabled.")
    find_package(OpenMP)

    if (OPENMP_FOUND)

        # For CMake < 3.9, we need to make the target ourselves
        if(NOT TARGET OpenMP::OpenMP_CXX)
            find_package(Threads REQUIRED)
            add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
            set_property(TARGET OpenMP::OpenMP_CXX
                    PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
            # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
            set_property(TARGET OpenMP::OpenMP_CXX
                    PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

        endif()
    else()
        message(FATAL_ERROR "Could not find OpenMP!")
    endif()
else()
    message(STATUS "OpenMP disabled.")
endif()
