message(STATUS "highway")

include(FetchContent)

FetchContent_Declare(
    autopas_highway
    URL ${AUTOPAS_SOURCE_DIR}/libs/highway-1.3.0.zip
)

FetchContent_GetProperties(autopas_highway)
if (NOT autopas_highway_POPULATED)
    FetchContent_Populate(autopas_highway)
    set(HWY_ENABLE_TESTS OFF CACHE BOOL "Disable tests")
    add_subdirectory(${autopas_highway_SOURCE_DIR} ${autopas_highway_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()