include(FetchContent)

FetchContent_Declare(
    autopas_highway
    SOURCE_DIR ${AUTOPAS_SOURCE_DIR}/libs/highway
)

# We must force this to ensure that Highway does not create a gtest target in conflict with our own.
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Disable Highway tests" FORCE)
FetchContent_MakeAvailable(autopas_highway)