message(STATUS "Using bundled Eigen3")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
    eigen3
    URL
        # eigen-master:
        # https://bitbucket.org/eigen/eigen/get/default.zip
        # eigen-3.3.7:
        ${CMAKE_SOURCE_DIR}/libs/eigen-eigen-323c052e1731.zip
    URL_HASH MD5=0d9c8496922d5c07609b9f3585f00e49
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen-3
    # since we only unpack a header lib src == include
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen-3/include
    # Disable build & install steps
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(eigen3 source_dir)

target_include_directories(autopas SYSTEM PUBLIC ${source_dir})
