message(STATUS "highway")

# Enable FetchContent CMake module
include(FetchContent)

# Build Eigen3 and make the cmake targets available
FetchContent_Declare(
        highway
        URL
        # highway
        ${AUTOPAS_SOURCE_DIR}/libs/highway
)

# Check if population has already been performed
FetchContent_GetProperties(highway)
if (NOT highway_POPULATED) # must be lowercase "eigen3" Fetch the content using previously declared
    # details
    FetchContent_Populate(highway)

    add_library(
            highway
            OBJECT
            IMPORTED GLOBAL
    )

    target_include_directories(highway SYSTEM INTERFACE "${highway_SOURCE_DIR}/highway-master")

endif ()


