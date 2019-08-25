message(STATUS "xcerces-c - using bundled version")

# Enable ExternalProject CMake module
include(ExternalProject)

# Extract eigen3
ExternalProject_Add(
    xerces
    URL
        # xerces:
        #https://www-eu.apache.org/dist//xerces/c/3/sources/xerces-c-3.2.2.zip
        # xerces-c-3.2.2:
        ${CMAKE_SOURCE_DIR}/libs/xerces-c-3.2.2.zip
#        URL_HASH= 4d6936efedad787ab1719a9dcab273e7
    UPDATE_COMMAND ""
    UPDATE_DISCONNECTED true
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/xerces
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/xerces/src
    CONFIGURE_COMMAND
    BUILD_COMMAND
    INSTALL_COMMAND
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(xerces source_dir)

add_dependencies(md-flexible xerces)

target_include_directories(md-flexible SYSTEM PUBLIC ${source_dir})
