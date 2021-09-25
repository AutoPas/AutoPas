# cmake module for adding ALL

# Enable ExternalProject CMake module
include(FetchContent)
FetchContent_Declare(
  all
  URL ${CMAKE_CURRENT_SOURCE_DIR}/ALL_20210925.zip
  URL_HASH MD5=75384d1773e28abae1b81b5b9cefe991
)

# Get ALL source and binary directories from CMake project
FetchContent_GetProperties(all)

if (NOT all_POPULATED)
  FetchContent_MakeAvailable(all)
endif ()
