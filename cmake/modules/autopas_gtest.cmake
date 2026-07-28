message(STATUS "gtest - using bundled version")
find_package(Threads REQUIRED)

option(INSTALL_GTEST "" OFF)

# hide options from ccmake
mark_as_advanced(
        BUILD_GMOCK
        INSTALL_GTEST
        GTEST_HAS_ABSL  # Only relevant when also using the Abseil library, which we don't.
)

# Prevent overriding the parent project's compiler/linker settings on Windows.
# => Compiles gtest with correct mt(d)/md(d)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# EXCLUDE_FROM_ALL, such that gtest is not built by default (target `all`).
add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/googletest ${CMAKE_CURRENT_BINARY_DIR}/googletest EXCLUDE_FROM_ALL)
