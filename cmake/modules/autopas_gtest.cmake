# check whether gtest is installed

message(STATUS "Using bundled gtest")
find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Download and install GoogleTest
ExternalProject_Add(
    gtest
    URL
        # https://github.com/google/googletest/archive/master.zip
        ${CMAKE_SOURCE_DIR}/libs/googletest-master.zip
    URL_HASH MD5=9ead2b6ec99010eb7ec77fdaf6d9ded9
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/gtest/src/gtest-build/lib/libgtest.a
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/gtest/src/gtest-build/lib/libgmock.a
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    # Disable install step
    INSTALL_COMMAND ""
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest STATIC IMPORTED GLOBAL)

add_dependencies(libgtest gtest)

# Set libgtest properties
set_target_properties(
    libgtest
    PROPERTIES
        "IMPORTED_LOCATION"
        "${binary_dir}/lib/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES"
        "${CMAKE_THREAD_LIBS_INIT}"
)

# Create a libgmock target to be used as a dependency by test programs
add_library(libgmock STATIC IMPORTED GLOBAL)

add_dependencies(libgmock gtest)

# Set libgmock properties
set_target_properties(
    libgmock
    PROPERTIES
        "IMPORTED_LOCATION"
        "${binary_dir}/lib/libgmock.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES"
        "${CMAKE_THREAD_LIBS_INIT}"
)

# I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
include_directories(SYSTEM "${source_dir}/googletest/include" "${source_dir}/googlemock/include")
