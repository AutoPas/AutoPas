# Enable FetchContent CMake module

include(FetchContent)

# Build Highway and make the cmake targets available
FetchContent_Declare(
    highway
    URL ${AUTOPAS_SOURCE_DIR}/libs/highway-master.zip
)

option(AUTOPAS_BUILD_HIGHWAY "" ON)

FetchContent_GetProperties(highway)
if (NOT highway_POPULATED)

    FetchContent_Populate(highway)
    message(STATUS "Highway")
    add_subdirectory(${highway_SOURCE_DIR} ${highway_BINARY_DIR})

endif()