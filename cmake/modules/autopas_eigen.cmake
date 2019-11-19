# capital E actually required...
find_package(Eigen3 QUIET)
# actually I don't know our minimal supported version but this is the one I tested.
if (Eigen3_FOUND AND "${Eigen3_VERSION}" VERSION_GREATER_EQUAL 3.3.7)
    message(STATUS "Eigen3 - using installed system version ${Eigen3_VERSION}")
    target_link_libraries(autopas PUBLIC Eigen3::Eigen)
else ()
    # system version not found -> install bundled version
    message(STATUS "Eigen3 - not found or version older than 3.3.7")
    message(
        STATUS
            "Eigen3 - if you want to use your version point the cmake variable Eigen3_DIR to the directory containing Eigen3Config.cmake in order to pass hints to find_package"
    )
    message(STATUS "Eigen3 - using bundled version 3.3.90 (commit 66be6c7)")

    # Enable ExternalProject CMake module
    include(ExternalProject)

    # Extract Eigen3
    ExternalProject_Add(
        Eigen3
        URL
            # eigen-master:
            # https://bitbucket.org/eigen/eigen/get/default.zip
            # eigen-3.3.90:
            ${CMAKE_SOURCE_DIR}/libs/eigen-eigen-66be6c76fc01.zip
        URL_HASH MD5=faaf36185ad92b039f7b3f641340dc28
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen-3
        # since we only unpack a header lib src == include
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen-3/include
        # Disable build & install steps
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )

    # Get GTest source and binary directories from CMake project
    ExternalProject_Get_Property(Eigen3 source_dir)

    add_dependencies(autopas Eigen3)

    target_include_directories(autopas SYSTEM PUBLIC ${source_dir})
endif ()
