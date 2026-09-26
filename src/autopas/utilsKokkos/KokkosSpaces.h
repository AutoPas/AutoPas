/**
 * @file KokkosSpaces.h
 * @author Luis Gall
 * @date 22. Sep 2026
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "Kokkos_Core.hpp"

namespace autopas::utilsKokkos {

using DeviceExecSpace = Kokkos::DefaultExecutionSpace;

using DeviceMemSpace = DeviceExecSpace::memory_space;

using HostExecSpace = Kokkos::HostSpace::execution_space;

using HostMemSpace = HostExecSpace::memory_space;

}

#endif
