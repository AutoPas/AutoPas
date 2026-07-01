option(AUTOPAS_ENABLE_KOKKOS "Enables the Kokkos containers and traversals" OFF)

if (NOT ${AUTOPAS_ENABLE_KOKKOS})
    return()
endif ()

set(Kokkos_VERSION 5.1.1)

find_package(Kokkos ${Kokkos_VERSION} CONFIG QUIET)

if (Kokkos_FOUND)
    message(STATUS "Found local Kokkos ${Kokkos_VERSION} Installation")
    return()
endif ()

message(STATUS "Using Kokkos from GitHub Release ${Kokkos_VERSION}")

include(FetchContent)

FetchContent_Declare(
        Kokkos
        URL
        https://github.com/kokkos/kokkos/archive/refs/tags/${Kokkos_VERSION}.tar.gz
)
FetchContent_MakeAvailable(Kokkos)