project(AutoPas)

cmake_minimum_required(VERSION 3.3)


#CMAKE policies, etc.
CMAKE_POLICY(SET CMP0057 NEW)
set(CMAKE_CXX_STANDARD 11)
# cmake < 3.6 ignores "set(CMAKE_CXX_STANDARD 11)" when using icpc
if (CMAKE_VERSION VERSION_LESS "3.6" AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()


set(CMAKE_INCLUDE_CURRENT_DIR ON)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

ENABLE_TESTING()


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)

include(other-compileroptions)
include(vectorization)
include(clang-sanitizers)
include(version.cmake)
include(clang-tidy)
include(clang-format)
include(doxygen)
include(openmp)
include(CodeCoverageAutopas)


#autopas
ADD_SUBDIRECTORY(src)

#libraries
ADD_SUBDIRECTORY(libs)

#tests and examples
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(tests)