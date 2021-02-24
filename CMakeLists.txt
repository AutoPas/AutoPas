cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(AutoPas)

# CMAKE policies, etc.
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
cmake_policy(SET CMP0079 NEW)

# Always use AUTOPAS_[SOURCE/BINARY]_DIR instead of CMAKE_[SOURCE/BINARY]_DIR or use
# CMAKE_CURRENT_[SOURCE/BINARY]_DIR directly for relative paths to allow including AutoPas via
# FetchContent as a subproject otherwise paths break because CMAKE_[SOURCE/BINARY]_DIR is defined by
# the parent project
set(AUTOPAS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(AUTOPAS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# starting with cmake 3.17 this is set using target_compile_features.
if (CMAKE_VERSION VERSION_LESS 3.17)
    message(
        STATUS
            "Setting CMAKE_CXX_STANDARD and CMAKE_CUDA_STANDARD globally because cmake < 3.17 is used."
    )
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CUDA_STANDARD 14)
endif ()

if (NOT CYGWIN)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif ()

set(SUPPORTED_COMPILERS "GNU;Intel;Clang")

if (NOT CMAKE_CXX_COMPILER_ID IN_LIST SUPPORTED_COMPILERS)
    message(WARNING "Compiler " ${CMAKE_CXX_COMPILER_ID} " not supported!")
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
else ()
    message(STATUS "Cuda disabled.")
endif ()

enable_testing()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

include(ccache)
include(clang-format)
include(clang-sanitizers)
include(clang-tidy)
include(cmake-format)
include(format)
include(CodeCoverageAutopas)
include(OpenMP)
include(mpi)

include(version.cmake)

# libraries
add_subdirectory(libs)


include(FetchContent)
FetchContent_Declare(
        kokkos
        GIT_REPOSITORY "https://github.com/kokkos/kokkos.git"
        GIT_TAG "3.3.01"
)
FetchContent_MakeAvailable(kokkos)


# autopas
add_subdirectory(src)
# autopasTools
add_subdirectory(tools)

# modules needing targets:
include(coloring)
include(logging)
include(other-compileroptions)
include(vectorization)
include(cuda)

# tests and examples
add_subdirectory(examples)
add_subdirectory(tests)

include(doxygen) # does not really need a target but the variable AUTOPAS_BUILD_EXAMPLES, which is defined in examples.

# write version information
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    AutoPasConfigVersion.cmake
    VERSION ${AutoPas_VERSION}
    COMPATIBILITY ExactVersion
    # for now ExactVersion is used, as the interface might change.
)
