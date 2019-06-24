message(STATUS "Using bundled Eigen3")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
        eigen3
        URL
            # https://bitbucket.org/eigen/eigen/get/default.zip
            ${CMAKE_SOURCE_DIR}/libs/eigen-eigen-323c052e1731.zip
        URL_HASH MD5=0d9c8496922d5c07609b9f3585f00e49
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen-3.3.7
        # Disable build & install steps
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/eigen-3.3.7/src/eigen3)