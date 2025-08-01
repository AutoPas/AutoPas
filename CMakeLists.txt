cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(AutoPas)

# CMAKE policies, etc.
cmake_policy(SET CMP0054 NEW) # from version 3.1
cmake_policy(SET CMP0057 NEW) # from version 3.3
cmake_policy(SET CMP0079 NEW) # from version 3.13
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.24)
    cmake_policy(SET CMP0135 NEW)
endif ()

# Always use AUTOPAS_[SOURCE/BINARY]_DIR instead of CMAKE_[SOURCE/BINARY]_DIR or use
# CMAKE_CURRENT_[SOURCE/BINARY]_DIR directly for relative paths to allow including AutoPas via
# FetchContent as a subproject otherwise paths break because CMAKE_[SOURCE/BINARY]_DIR is defined by
# the parent project
set(AUTOPAS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(AUTOPAS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# append custom modules to the module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)
# load custom conditional options
include(optionConditional)

# Determine if AutoPas is built as the primary project or not.
option_conditional(AUTOPAS_STANDALONE_BUILD
        "Add options and targets which are only relevant for developers of AutoPas"
        "AUTOPAS_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR"
)

# starting with cmake 3.17 this is set using target_compile_features.
if (CMAKE_VERSION VERSION_LESS 3.17)
    message(
            STATUS
            "Setting CMAKE_CXX_STANDARD globally because cmake < 3.17 is used."
    )
    set(CMAKE_CXX_STANDARD 20)
endif ()

if (NOT CYGWIN)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif ()

set(SUPPORTED_COMPILERS "GNU;IntelLLVM;Clang;AppleClang")

if (NOT CMAKE_CXX_COMPILER_ID IN_LIST SUPPORTED_COMPILERS)
    message(WARNING "Compiler " ${CMAKE_CXX_COMPILER_ID} " not supported!")
endif ()

if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin" AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    # As of macOS 15.4, _Alignof doesn't seem to be specified anymore in the macOS system headers
    # leading to compile failure when using GCC instead of Apple Clang. This is a work-around.
    add_definitions(-D_Alignof=alignof)
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
endif ()

enable_testing()

option_conditional(AUTOPAS_FORMATTING_TARGETS
        "Adds targets for code formatting in AutoPas. Highly useful when developing AutoPas but might cause conflicts when only using it."
        AUTOPAS_STANDALONE_BUILD
)

if (AUTOPAS_FORMATTING_TARGETS)
    include(autopas_clang-format)
    include(autopas_cmake-format)
    include(autopas_format)
endif ()

include(autopas_ccache)
include(autopas_clang-sanitizers)
include(autopas_clang-tidy)
include(autopas_OpenMP)
include(autopas_mpi)

include(version.cmake)

# libraries
add_subdirectory(libs)

include(autopas_antlr4cpp)

# autopas
add_subdirectory(src)
# autopasTools
add_subdirectory(tools)

# modules needing targets:
include(autopas_coloring)
include(autopas_logging)
include(autopas_other-compileroptions)
include(autopas_vectorization)

# tests and examples
add_subdirectory(tests)
# application libraries
add_subdirectory(applicationLibrary)
# might also include tests, hence might need flags from tests subdirectory
add_subdirectory(examples)

include(autopas_doxygen) # does not really need a target but the variable AUTOPAS_BUILD_EXAMPLES, which is defined in examples.

# write version information
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        AutoPasConfigVersion.cmake
        VERSION ${AutoPas_VERSION}
        COMPATIBILITY ExactVersion
        # for now ExactVersion is used, as the interface might change.
)

# We have to include coverage target after the CMakeFiles for the tests
include(autopas_CodeCoverage)