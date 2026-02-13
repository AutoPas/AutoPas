include(FetchContent)

FetchContent_Declare(
    autopas_highway
    URL ${AUTOPAS_SOURCE_DIR}/libs/highway-1.3.0.zip
)

# We must force this to ensure that Highway does not create a gtest target in conflict with our own.
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Disable Highway tests" FORCE)
FetchContent_MakeAvailable(autopas_highway)