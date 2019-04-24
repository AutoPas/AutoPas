option(KOKKOS_ENABLED "Activate Kokkos Options." OFF)
#option(KOKKOS_TARGET "Please choose your device." OFF)


set_property(CACHE KOKKOS_DEVICE PROPERTY STRINGS "Serial;OpenMP;Cuda")

# set Release as the default build type if it is not yet set.
if(NOT KOKKOS_DEVICE)
    set(KOKKOS_DEVICE "Serial" CACHE STRING
            "Choose the appropiate kokkos_device, options are: Serial OpenMP Cuda." FORCE)
endif(NOT KOKKOS_DEVICE)
#look for kokkos

#include kokkos

if(KOKKOS_ENABLED)

    message("Kokkos is enabled")

    #set path to kokkos files
    set(KOKKOS_INSTALL_PATH "$ENV{HOME}/kokkos")
    set(KOKKOS_INSTALL_PATH_BIN "${HOME}/kokkos/build")

    message("Kokkos_Device: " ${KOKKOS_DEVICE})

    #add_subdirectory(${KOKKOS_INSTALL_PATH} ${KOKKOS_INSTALL_PATH_BIN})
    #include_directories(${Kokkos_INCLUDE_DIRS_RET})

    if(${KOKKOS_DEVICE} STREQUAL "OpenMP")
        message("Kokkos uses OpenMP")
        set(KOKKOS_ENABLE_OPENMP yes)
        set(OMP_PROC_BIND=spread)
        set(OMP_PLACES threads)
    elseif(${KOKKOS_DEVICE}STREQUAL Cuda)
        message("TODO Cuda is not yet implemented, Serial set")
    else()
        message("Kokkos uses Serial (no parallelization)")
    endif()
else()
    message("Kokkos is disabled")
endif()