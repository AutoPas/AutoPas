option(AUTOPAS_KOKKOS "Activate Kokkos in Autopas." OFF)

#set_property(CACHE KOKKOS_DEVICE PROPERTY STRINGS "Serial;OpenMP;Cuda")

# set Release as the default build type if it is not yet set.
#if(NOT KOKKOS_DEVICE)
#    set(KOKKOS_DEVICE "Serial" CACHE STRING
#            "Choose the appropiate kokkos_device, options are: Serial OpenMP Cuda." FORCE)
#endif(NOT KOKKOS_DEVICE)

if(AUTOPAS_KOKKOS)
    message(STATUS "Kokkos enabled.")

    #set path to kokkos files
    set(KOKKOS_INSTALL_PATH "$ENV{HOME}/kokkos" CACHE STRING "Path to kokkos project root directory")
    #check, whether directory to kokkos exists, no check on functionality of the content
    if(EXISTS "${KOKKOS_INSTALL_PATH}" AND IS_DIRECTORY "${KOKKOS_INSTALL_PATH}")
        message(STATUS "Kokkos Directory found")
    else()
        message(FATAL_ERROR "Could not find Kokkos! Please set KOKKOS_INSTALL_PATH to the kokkos directory!")
    endif()
    set(CMAKE_CXX_EXTENSIONS OFF)

    # build kokkos
    add_subdirectory(${KOKKOS_INSTALL_PATH} ${PROJECT_BINARY_DIR}/kokkos)

    # if user indicates to use OpenMP set all variables
    # @FIXME this deletes the docstrings of the updated options
    if(${KOKKOS_ENABLE_OPENMP} OR ${OPENMP})
        set(KOKKOS_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
        set(OPENMP ON CACHE BOOL "" FORCE)
        message(STATUS "Kokkos uses OpenMP.")
    # if user indicates to use CUDA set all variables
    elseif(${KOKKOS_ENABLE_CUDA}) # OR ${ENABLE_CUDA})
        set(ENABLE_CUDA ON CACHE BOOL FORCE)
        set(KOKKOS_ENABLE_CUDA ON CACHE BOOL FORCE)
        set(Kokkos_ENABLE_Cuda_Lambda ON CACHE BOOL FORCE)
        message(STATUS "Kokkos uses Cuda with KOKKOS_ARCH: " ${KOKKOS_ARCH})
    else()
        message(WARNING "Kokkos is in serial mode (no parallelization).")
    endif()

endif()
