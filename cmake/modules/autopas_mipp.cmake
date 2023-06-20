# Enable FetchContent CMake module
include(FetchContent)

# Build Eigen3 and make the cmake targets available
FetchContent_Declare(
    mipp
    URL
        ${AUTOPAS_SOURCE_DIR}/libs/MIPP-master.zip
    #URL_HASH MD5=994092410ba29875184f7725e0371596
)

# Check if population has already been performed
FetchContent_GetProperties(mipp)
if (NOT mipp_POPULATED) # must be lowercase "eigen3" Fetch the content using previously declared
                          # details
    FetchContent_Populate(mipp)


    message(STATUS "MIPP")
    add_library(
            mipp
        OBJECT # this is a header only lib therefore object type is needed
        IMPORTED GLOBAL
    )
    target_include_directories(mipp SYSTEM INTERFACE "${mipp_SOURCE_DIR}")

    # add_subdirectory(${eigen3_SOURCE_DIR} ${eigen3_BINARY_DIR})
endif ()
