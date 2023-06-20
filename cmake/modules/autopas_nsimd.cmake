message(STATUS "nsimd")

# Enable FetchContent CMake module
include(FetchContent)

# Build Eigen3 and make the cmake targets available
FetchContent_Declare(
        nsimd
        URL
        # nsimd
        ${AUTOPAS_SOURCE_DIR}/libs/nsimd-master.zip
)

# Check if population has already been performed
FetchContent_GetProperties(nsimd)
if (NOT nsimd_POPULATED) # must be lowercase "eigen3" Fetch the content using previously declared
    # details
    FetchContent_Populate(nsimd)

    # Do not add_subdirectory, else we would run configure, build and install Just define a library
    # from the sources
    add_library(
            nsimd
            OBJECT # this is a header only lib therefore object type is needed
            IMPORTED GLOBAL
    )

    target_include_directories(nsimd SYSTEM INTERFACE "${nsimd_SOURCE_DIR}/include")
endif ()