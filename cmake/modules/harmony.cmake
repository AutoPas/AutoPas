IF (MSVC)
    message(STATUS "MSVC detected, disabling active harmony support.")
    return()
ENDIF()

message(STATUS "harmony - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# harmony is a make file project
find_program(MAKE_EXE NAMES gmake nmake make)

# Prefer using ExternalProject for harmony over FetchContent. FetchContent is only really useful for
# other cmake projects as it allows the export of their cmake targets. As harmony is a Makefile
# project, FetchContent is not useful here.

# Extract and build harmony
ExternalProject_Add(
    harmony_bundled
    URL ${AUTOPAS_SOURCE_DIR}/libs/harmony.zip
    URL_HASH MD5=a8768c2886bdc2e44e3b6b7d4f94729c
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/harmony/include/lib/libharmony.a
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/harmony
    # since we only unpack a header lib src == include
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/harmony/include
    # tell cmake to run make inside the source folder and suppress all warnings
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    LOG_BUILD ON
    LOG_INSTALL ON
    BUILD_COMMAND ${MAKE_EXE}
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(harmony_bundled source_dir)

add_library(
    harmony
    STATIC
    IMPORTED
    GLOBAL
)

add_dependencies(harmony harmony_bundled)

set_target_properties(
    harmony
    PROPERTIES "IMPORTED_LOCATION" "${CMAKE_CURRENT_BINARY_DIR}/harmony/include/lib/libharmony.a"
)

# create directory otherwise cmake will complain during generate step bc it would only be created by
# make
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/harmony/include/include")

target_include_directories(
    harmony SYSTEM
    INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/harmony/include/include"
)

# Set macro needed to set environment variable for ActiveHarmony
target_compile_definitions(
    harmony
        INTERFACE
        HARMONY_HOME="HARMONY_HOME=${CMAKE_CURRENT_BINARY_DIR}/harmony/include"
        AUTOPAS_ENABLE_ACTIVE_HARMONY
)
