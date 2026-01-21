option(AUTOPAS_ENABLE_KOKKOS "Enables the GPU backend using Kokkos" OFF)

set(AUTOPAS_KOKKOS_BACKEND "CUDA" CACHE STRING "Selects Kokkos device backend (CUDA, HIP, or SYCL, OPENMP)" )

set_property(CACHE AUTOPAS_KOKKOS_BACKEND PROPERTY STRINGS CUDA HIP SYCL OPENMP)


if (NOT ${AUTOPAS_ENABLE_KOKKOS})
    return()
endif ()


set(Kokkos_ENABLE_OPENMP ON)
set(Kokkos_ENABLE_SERIAL ON)

if (${AUTOPAS_KOKKOS_BACKEND} STREQUAL "CUDA")
    set(Kokkos_ENABLE_CUDA ON)
    #set(Kokkos_ARCH_AMPERE86 ON)
    set(Kokkos_ARCH_PASCAL61 ON)
    set(Kokkos_ENABLE_CUDA_CONSTEXPR ON)

elseif (${AUTOPAS_KOKKOS_BACKEND} STREQUAL "HIP")
    set(Kokkos_ENABLE_HIP ON)

elseif (${AUTOPAS_KOKKOS_BACKEND} STREQUAL "SYCL")
    set(Kokkos_ENABLE_SYCL ON)
    set(Kokkos_ARCH_INTEL_GEN12LP ON)
elseif (${AUTOPAS_KOKKOS_BACKEND} STREQUAL "OPENMP")
else ()
    message(FATAL_ERROR "Unsupported Kokkos backend selected: ${AUTOPAS_KOKKOS_BACKEND}. Supported are CUDA, HIP, SYCL, OPENMP.")
endif ()

message(STATUS "Setting up Kokkos")
set(Kokkos_VERSION 4.7.01)

#find_package(Kokkos ${Kokkos_VERSION} CONFIG QUIET)

if (${Kokkos_FOUND})
    message(STATUS "Found existing Kokkos libraries: ${Kokkos_DIR}")
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
