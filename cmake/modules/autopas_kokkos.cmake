# oriented on branch from JoHaHo at commit id 9c27df2

set(Kokkos_ENABLE_OPENMP ON)
set(Kokkos_ENABLE_CUDA ON)
set(Kokkos_ENABLE_CUDA_CONSTEXPR ON)

set(SPDLOG_USE_STD_FORMAT ON)
set(AUTOPAS_MIN_LOG_LEVEL OFF)
set(AUTOPAS_OPENMP OFF)

set(Kokkos_ARCH_PASCAL61 ON)

set(Kokkos_VERSION 4.7.01)

message("Fetching Kokkos from Github")

include(FetchContent)

FetchContent_Declare(
        Kokkos
        URL https://github.com/kokkos/kokkos/archive/refs/tags/${Kokkos_VERSION}.tar.gz
)
FetchContent_MakeAvailable(Kokkos)