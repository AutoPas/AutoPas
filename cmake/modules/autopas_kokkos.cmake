option(AUTOPAS_ENABLE_KOKKOS "Enables the Kokkos containers and traversals" OFF)

if (NOT AUTOPAS_ENABLE_KOKKOS)
    return()
endif ()

message(STATUS "Setting up Kokkos")
set(Kokkos_VERSION 5.1.1)

# Detect which accelerator compilers are available. This has to happen before the Kokkos backend is
# selected, as the backend selection below inspects the enabled languages.
include(CheckLanguage)

check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit QUIET)
    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
        # "native" detects the architecture of the GPU of the build machine.
        set(CMAKE_CUDA_ARCHITECTURES "native" CACHE STRING "CUDA architectures to compile for")
    else ()
        message(WARNING "CMake < 3.24 cannot detect CUDA architectures natively. Set CMAKE_CUDA_ARCHITECTURES manually.")
    endif ()
    if (CUDAToolkit_FOUND)
        message(STATUS "CUDA language enabled and CUDAToolkit found")
    else ()
        message(WARNING "CUDA compiler available but CUDAToolkit not found")
    endif ()
else ()
    message(STATUS "CUDA compiler not available with current toolchain")
endif ()

check_language(HIP)
if (CMAKE_HIP_COMPILER)
    enable_language(HIP)
    find_package(HIP QUIET)
    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
        # "native" detects the architecture of the GPU of the build machine.
        set(CMAKE_HIP_ARCHITECTURES "native" CACHE STRING "HIP architectures to compile for")
    else ()
        message(WARNING "CMake < 3.24 cannot detect HIP architectures natively. Set CMAKE_HIP_ARCHITECTURES manually.")
    endif ()
    if (HIP_FOUND)
        message(STATUS "HIP language enabled and HIP package found")
    else ()
        message(WARNING "HIP compiler available but HIP package not found")
    endif ()
else ()
    message(STATUS "HIP language not available with current toolchain")
endif ()

# Select a Kokkos backend according to the enabled languages, unless the user chose one explicitly.
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if ("CUDA" IN_LIST languages)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable Kokkos CUDA backend")
elseif ("HIP" IN_LIST languages)
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable Kokkos HIP backend")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(Kokkos_ENABLE_SYCL ON CACHE BOOL "Enable Kokkos SYCL backend")
else ()
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable Kokkos OpenMP backend")
endif ()

find_package(Kokkos ${Kokkos_VERSION} CONFIG QUIET)

if (Kokkos_FOUND)
    message(STATUS "Found local Kokkos ${Kokkos_VERSION} installation")
    return()
endif ()

message(STATUS "Using Kokkos from GitHub Release ${Kokkos_VERSION}")

include(FetchContent)

# Always build for the machine on which is being compiled, unless specified otherwise.
set(Kokkos_ARCH_NATIVE ON CACHE BOOL "Build Kokkos for the machine on which it is being compiled")

FetchContent_Declare(
        Kokkos
        URL https://github.com/kokkos/kokkos/archive/refs/tags/${Kokkos_VERSION}.tar.gz
        URL_HASH SHA256=77cbde0066f5ea9343d35be452826b6b226ceacb385239c28cc9688baf471cc0
)
FetchContent_MakeAvailable(Kokkos)
