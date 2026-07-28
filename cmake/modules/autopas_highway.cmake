# We must force this to ensure that Highway does not create a gtest target in conflict with our own.
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Disable Highway tests" FORCE)

add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/highway ${CMAKE_CURRENT_BINARY_DIR}/highway EXCLUDE_FROM_ALL)