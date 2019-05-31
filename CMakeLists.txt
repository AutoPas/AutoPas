project(AutoPas)

cmake_minimum_required(VERSION 3.10)


#CMAKE policies, etc.
CMAKE_POLICY(SET CMP0054 NEW)
CMAKE_POLICY(SET CMP0057 NEW)
set(CMAKE_CXX_STANDARD 14)

set(SUPPORTED_COMPILERS "GNU;Intel;Clang")

if (NOT CMAKE_CXX_COMPILER_ID IN_LIST SUPPORTED_COMPILERS)
    message(WARNING "Compiler " ${CMAKE_CXX_COMPILER_ID} " not supported!")
endif ()

# Check if the compiler command is archer
# CMAKE_CXX_COMPILER contains the whole path. "archer" must be after the last "/"
if (CMAKE_CXX_COMPILER MATCHES ".*archer[^/]*")
    message(STATUS "Archer detected!")
    set(ARCHER true)
endif ()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# CUDA options
option(ENABLE_CUDA "Activates Cuda Components" OFF)
if(ENABLE_CUDA)
    message(STATUS "Cuda enabled.")
    enable_language(CUDA)
    set(CUDA_SEPARABLE_COMPILATION ON)
else()
    message(STATUS "Cuda disabled.")
endif()

ENABLE_TESTING()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)

include(ccache)
include(clang-format)
include(clang-sanitizers)
include(clang-tidy)
include(CodeCoverageAutopas)
include(doxygen)
include(OpenMP)
include(coloring)

include(version.cmake)

#autopas
ADD_SUBDIRECTORY(src)

#libraries
ADD_SUBDIRECTORY(libs)

#tests and examples
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(tests)

# modules needing targets:
include(logging)
include(other-compileroptions)
include(vectorization)
include(cuda)
