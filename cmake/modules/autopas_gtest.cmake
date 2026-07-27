message(STATUS "gtest - using bundled version")
find_package(Threads REQUIRED)

# Enable FetchContent CMake module
include(FetchContent)

# Build GoogleTest and make the cmake targets available
FetchContent_Declare(
        autopas_gtest
        SOURCE_DIR ${AUTOPAS_SOURCE_DIR}/libs/googletest
)

option(INSTALL_GTEST "" OFF)

# hide options from ccmake
mark_as_advanced(BUILD_GMOCK INSTALL_GTEST)

# Prevent overriding the parent project's compiler/linker settings on Windows.
# => Compiles gtest with correct mt(d)/md(d)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(autopas_gtest)
