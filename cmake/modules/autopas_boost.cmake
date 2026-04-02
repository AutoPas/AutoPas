option(Boost_ForceBundled "Do not look for an installed version, always use bundled." ON)


message(STATUS "Fast particle buffer with linear regression estimator: enabling Boost")
if (NOT ${Boost_ForceBundled})
    find_package(Boost QUIET)
    if (Boost_FOUND)
        message(STATUS "Boost - using installed system version ${Boost_VERSION}")
        return()
    else ()
        message(STATUS "Boost - installed system version not found using bundled one.")
    endif ()
endif ()

include(FetchContent)
FetchContent_Declare(
        boost
        URL ${CMAKE_SOURCE_DIR}/libs/boost-1.90.0.zip
        URL_HASH MD5=dcf37834eb71f7a6e49482b226e246a1
)

FetchContent_MakeAvailable(boost)

add_library(boost INTERFACE)

target_include_directories(boost SYSTEM INTERFACE
        ${boost_SOURCE_DIR}/include
)

