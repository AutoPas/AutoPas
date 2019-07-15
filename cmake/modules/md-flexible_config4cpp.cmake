message(STATUS "config4cpp - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
        config4cpp
        URL
        # config4cpp:
        # http://www.config4star.org/download/config4cpp.zip
        # config4cpp:
        ${CMAKE_SOURCE_DIR}/libs/config4cpp.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/config4cpp
        # since we only unpack a header lib src == include
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/config4cpp/include
        # Disable build & install steps
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(config4cpp source_dir)

add_dependencies(md-flexible config4cpp)

target_include_directories(autopas SYSTEM PUBLIC ${source_dir})
