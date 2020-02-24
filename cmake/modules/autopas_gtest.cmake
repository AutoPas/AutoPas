message(STATUS "gtest - using bundled version")
find_package(Threads REQUIRED)

# Enable FetchContent CMake module
include(FetchContent)

# Build GoogleTest and make the cmake targets available
FetchContent_Declare(
    gtest
    URL ${AUTOPAS_SOURCE_DIR}/libs/googletest-1.10.0.zip
    URL_HASH MD5=82358affdd7ab94854c8ee73a180fc53
)

option(INSTALL_GTEST "" OFF)

# hide options from ccmake
mark_as_advanced(BUILD_GMOCK INSTALL_GTEST)

FetchContent_GetProperties(gtest)
if (NOT gtest_POPULATED)
    FetchContent_Populate(gtest)
    add_subdirectory(${gtest_SOURCE_DIR} ${gtest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
