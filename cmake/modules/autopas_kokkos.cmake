include(FetchContent)

FetchContent_Declare(
    Kokkos
    URL https://github.com/kokkos/kokkos/archive/refs/tags/4.5.01.zip
    # URL_HASH TODO: <find in Github>
    )

FetchContent_MakeAvailable(Kokkos)