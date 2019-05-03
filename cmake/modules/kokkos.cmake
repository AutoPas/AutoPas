option(KOKKOS_ENABLED "Activate Kokkos in Autopas." OFF)



#set_property(CACHE KOKKOS_DEVICE PROPERTY STRINGS "Serial;OpenMP;Cuda")

# set Release as the default build type if it is not yet set.
#if(NOT KOKKOS_DEVICE)
#    set(KOKKOS_DEVICE "Serial" CACHE STRING
#            "Choose the appropiate kokkos_device, options are: Serial OpenMP Cuda." FORCE)
#endif(NOT KOKKOS_DEVICE)



if(KOKKOS_ENABLED)

    message("Kokkos is enabled; KOKKOS_DEVICE: " ${KOKKOS_DEVICE})
    message("Compiler: " ${CMAKE_CXX_COMPILER})
    #set path to kokkos files
    set(KOKKOS_INSTALL_PATH "$ENV{HOME}/kokkos")
    set(CMAKE_CXX_EXTENSIONS OFF)

    #message("Kokkos_Device: " ${KOKKOS_DEVICE})

    #if(${KOKKOS_ENABLED} STREQUAL "ON")

    add_definitions(-DKOKKOS_ENABLED="on" )

    set(KOKKOS_DEVICE "Serial") #default: Serial

    add_subdirectory(${KOKKOS_INSTALL_PATH} ${PROJECT_BINARY_DIR}/kokkos)
    include_directories(${Kokkos_INCLUDE_DIRS_RET})

    if(${KOKKOS_DEVICE} STREQUAL OpenMP)
        message("Kokkos uses OpenMP")
        set(KOKKOS_ENABLE_OPENMP yes)
        set(OMP_PROC_BIND=spread)
        set(OMP_PLACES threads)
    elseif(${KOKKOS_DEVICE} STREQUAL Cuda)

        set(KOKKOS_ENABLE_CUDA ON)
        set(Kokkos_ENABLE_Cuda_Lambda ON)
        set(KOKKOS_ARCH "Kepler30")
        set(KOKKOS_DEVICES "Cuda")
        message("Kokkos uses Cuda with KOKKOS_ARCH: " ${KOKKOS_ARCH})
    else()
        message("Kokkos uses Serial (no parallelization)")
    endif()

else()
    message("Kokkos is disabled")
endif()