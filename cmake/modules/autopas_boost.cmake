if (AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER_LIN OR AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER)
    message(STATUS "Fast particle buffer with linear regression estimator: enabling Boost")
    #find_package(Boost REQUIRED)
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
else ()
    message(STATUS "Boost disabled.")
endif ()