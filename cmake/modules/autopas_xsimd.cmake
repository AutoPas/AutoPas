message(STATUS "xsimd")

# Enable FetchContent CMake module
include(FetchContent)

# Build Eigen3 and make the cmake targets available
FetchContent_Declare(
        xsimd
        URL
        # xsimd
        ${AUTOPAS_SOURCE_DIR}/libs/xsimd-master.zip
)

# Check if population has already been performed
FetchContent_GetProperties(xsimd)
if (NOT xsimd_POPULATED) # must be lowercase "eigen3" Fetch the content using previously declared
    # details
    FetchContent_Populate(xsimd)

    # Do not add_subdirectory, else we would run configure, build and install Just define a library
    # from the sources
    add_library(
            xsimd
            OBJECT # this is a header only lib therefore object type is needed
            IMPORTED GLOBAL
    )

    target_include_directories(xsimd SYSTEM INTERFACE "${xsimd_SOURCE_DIR}")
endif ()