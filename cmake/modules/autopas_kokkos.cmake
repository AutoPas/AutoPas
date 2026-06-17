option(AUTOPAS_ENABLE_KOKKOS "Enables the GPU backend using Kokkos" OFF)

if (NOT ${AUTOPAS_ENABLE_KOKKOS})
    return()
endif ()

set(Kokkos_VERSION 5.1.1)

find_package(Kokkos ${Kokkos_VERSION} CONFIG QUIET)

if (Kokkos_FOUND)
    message(STATUS "Found local Kokkos Installation")
else ()
    message(STATUS "Using Kokkos from GitHub Release ${Kokkos_VERSION}")
    include(FetchContent)

    # For the CPU Code always optimize for the machine being build on (use vectorization, etc.)
    set(Kokkos_ARCH_NATIVE ON CACHE BOOL "Always build for the machine on which AutoPas is being compiled" FORCE)

    FetchContent_Declare(
            Kokkos
            URL https://github.com/kokkos/kokkos/archive/refs/tags/${Kokkos_VERSION}.tar.gz
    )
    FetchContent_MakeAvailable(Kokkos)

    # Mark all CMake variables of the Kokkos project as advanced for this project expect the main backend selection
    get_cmake_property(_vars CACHE_VARIABLES)
    foreach (_var ${_vars})
        if (_var MATCHES "^Kokkos_" AND
                NOT _var MATCHES "^Kokkos_ENABLE_(SERIAL|OPENMP|THREADS|HPX|CUDA|HIP|SYCL|OPENMPTARGET|OPENACC)$")
            mark_as_advanced(${_var})
        endif ()
    endforeach ()

endif ()
