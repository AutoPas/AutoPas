option(AUTOPAS_ENABLE_KOKKOS "Enables the Kokkos containers and traversals" OFF)

if (NOT AUTOPAS_ENABLE_KOKKOS)
    return()
endif ()

message(STATUS "Setting up Kokkos")
set(Kokkos_VERSION 5.1.1)

get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if ("CUDA" IN_LIST languages)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable Kokkos CUDA backend" FORCE)
elseif ("HIP" IN_LIST languages)
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable Kokkos HIP backend" FORCE)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(Kokkos_ENABLE_SYCLE ON CACHE BOOL "Enable Kokkos Sycl backend" FORCE)
else ()
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable Kokkos OpenMP backend" FORCE)
endif ()

find_package(Kokkos ${Kokkos_VERSION} CONFIG QUIET)

if (Kokkos_FOUND)
    message(STATUS "Found local Kokkos ${Kokkos_VERSION} Installation")
    return()
endif ()

message(STATUS "Using Kokkos from GitHub Release ${Kokkos_VERSION}")

include(FetchContent)

set(Kokkos_ARCH_NATIVE ON CACHE STRING "Always build for the machine on which is being compiled" FORCE)

FetchContent_Declare(
        Kokkos
        URL
        https://github.com/kokkos/kokkos/archive/refs/tags/${Kokkos_VERSION}.tar.gz
)
FetchContent_MakeAvailable(Kokkos)