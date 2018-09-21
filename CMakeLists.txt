project(AutoPas)

cmake_minimum_required(VERSION 3.3)


#CMAKE policies, etc.
CMAKE_POLICY(SET CMP0057 NEW)
set(CMAKE_CXX_STANDARD 14)
# cmake < 3.6 ignores "set(CMAKE_CXX_STANDARD 11)" when using icpc
if (CMAKE_VERSION VERSION_LESS "3.6" AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif ()


set(CMAKE_INCLUDE_CURRENT_DIR ON)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# option for more verbose log messages.
option(VERBOSE_LOGGING "Print Filename and line in every log line." OFF)
if (VERBOSE_LOGGING)
    add_definitions(-DAUTOPAS_VERBOSE_LOG)
endif ()

ENABLE_TESTING()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)

include(other-compileroptions)
include(vectorization)
include(clang-sanitizers)
include(version.cmake)
include(clang-tidy)
include(clang-format)
include(doxygen)
include(OpenMP)
include(CodeCoverageAutopas)


#autopas
ADD_SUBDIRECTORY(src)

#libraries
ADD_SUBDIRECTORY(libs)

#tests and examples
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(tests)
