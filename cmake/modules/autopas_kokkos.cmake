include(FetchContent)

FetchContent_Declare(
    Kokkos
    URL https://github.com/kokkos/kokkos/releases/download/5.0.0/kokkos-5.0.0.tar.gz
    # URL_HASH TODO: <find in Github>
    )

FetchContent_MakeAvailable(Kokkos)