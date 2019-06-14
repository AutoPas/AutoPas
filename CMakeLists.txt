project(AutoPas)

cmake_minimum_required(VERSION 3.14)

# CMAKE policies, etc.
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
cmake_policy(SET CMP0079 NEW)

set(CMAKE_CXX_STANDARD 14)

set(SUPPORTED_COMPILERS "GNU;Intel;Clang")

if (NOT CMAKE_CXX_COMPILER_ID IN_LIST SUPPORTED_COMPILERS)
    message(WARNING "Compiler " ${CMAKE_CXX_COMPILER_ID} " not supported!")
endif ()

# ~~~
# Check if the compiler command is archer.
# CMAKE_CXX_COMPILER contains the whole path.
# "archer" must be after the last "/"
# ~~~
if (CMAKE_CXX_COMPILER MATCHES ".*archer[^/]*")
    message(STATUS "Archer detected!")
    set(ARCHER true)
endif ()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(
        CMAKE_BUILD_TYPE
        "Release"
        CACHE
            STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
            FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)

# CUDA options
option(AUTOPAS_ENABLE_CUDA "Activates Cuda Components" OFF)
if (AUTOPAS_ENABLE_CUDA)
    message(STATUS "Cuda enabled.")
    enable_language(CUDA)
    set(CUDA_SEPARABLE_COMPILATION ON)
else ()
    message(STATUS "Cuda disabled.")
endif ()

enable_testing()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)

include(clang-format)
include(clang-sanitizers)
include(clang-tidy)
include(cmake-format)
include(CodeCoverageAutopas)
include(doxygen)
include(OpenMP)

include(version.cmake)

# autopas
add_subdirectory(src)

# libraries
add_subdirectory(libs)

# tests and examples
add_subdirectory(examples)
add_subdirectory(tests)

# modules needing targets:
include(coloring)
include(logging)
include(other-compileroptions)
include(vectorization)
include(cuda)
