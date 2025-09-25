# cmake module for adding ALL

option(MD_FLEXIBLE_ENABLE_ALLLBL "Enable load balancing via ALL for MD-Flex" OFF)

# Enable ExternalProject CMake module
if (MD_FLEXIBLE_ENABLE_ALLLBL)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMD_FLEXIBLE_ENABLE_ALLLBL")

  include(FetchContent)
  FetchContent_Declare(
    allfetch
    URL ${AUTOPAS_SOURCE_DIR}/libs/ALL-0.9.3.zip
    URL_HASH MD5=aa30549b3decce8df1ee82e90b198b04
  )
  
  # Get ALL source and binary directories from CMake project
  FetchContent_GetProperties(allfetch)
  
  if (NOT allfetch)
    FetchContent_MakeAvailable(allfetch)
  endif ()
  
  set(ALL_LIB "ALL")
else()
  message(STATUS "ALL load balancing library support disabled")
  set(ALL_LIB "")
endif()
