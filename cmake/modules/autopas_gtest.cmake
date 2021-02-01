message(STATUS "gtest - using bundled version")
find_package(Threads REQUIRED)

# Enable FetchContent CMake module
include(FetchContent)

# Build GoogleTest and make the cmake targets available
FetchContent_Declare(
    autopas_gtest
    URL ${AUTOPAS_SOURCE_DIR}/libs/googletest-1.10.0.zip
    URL_HASH MD5=82358affdd7ab94854c8ee73a180fc53
)

option(INSTALL_GTEST "" OFF)

# hide options from ccmake
mark_as_advanced(BUILD_GMOCK INSTALL_GTEST)

# Prevent overriding the parent project's compiler/linker settings on Windows.
# => Compiles gtest with correct mt(d)/md(d)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_GetProperties(autopas_gtest)
if (NOT autopas_gtest_POPULATED)
    FetchContent_Populate(autopas_gtest)
    add_subdirectory(${autopas_gtest_SOURCE_DIR} ${autopas_gtest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()

target_compile_definitions(
        gmock
        PUBLIC
        $<$<BOOL:${MSVC}>:_ITERATOR_DEBUG_LEVEL=0>
)