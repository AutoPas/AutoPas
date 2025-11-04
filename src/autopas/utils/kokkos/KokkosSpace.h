#pragma once
#include <Kokkos_Core_fwd.hpp>

namespace autopas::utils::kokkos {

#ifdef KOKKOS_ENABLE_CUDA
using DeviceSpace = Kokkos::Cuda;
#elif KOKKOS_ENABLE_HIP
using DeviceSpace = Kokkos::Hip;
#endif

// using SharedSpace = Kokkos::SharedSpace;

#ifdef KOKKOS_ENABLE_OPENMP
using HostSpace = Kokkos::OpenMP;
#else
using HostSpace = Kokkos::Serial;
#endif

}  // namespace autopas::utils::kokkos
