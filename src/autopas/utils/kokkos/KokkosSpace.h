#pragma once
#include <Kokkos_Core_fwd.hpp>

namespace autopas::utils::kokkos {

#ifdef KOKKOS_ENABLE_CUDA
using DeviceSpace = Kokkos::Cuda;
#elif KOKKOS_ENABLE_HIP
using DeviceSpace = Kokkos::Hip;
#elif KOKKOS_ENABLE_SYCL
using DeviceSpace = Kokkos::Sycl;
#else
using DeviceSpace = Kokkos::OpenMP;
#endif

using SharedSpace = Kokkos::SharedSpace;

using HostSpace = Kokkos::DefaultHostExecutionSpace;

}  // namespace autopas::utils::kokkos
