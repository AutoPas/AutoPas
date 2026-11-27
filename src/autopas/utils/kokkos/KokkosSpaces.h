/**
 * @file KokkosSpaces.h
 * @author Luis Gall
 * @date 22. Sep 2026
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include "Kokkos_Core.hpp"

namespace autopas::utils::kokkos {

/**
 * Default Kokkos Device Execution Space. If there is no GPU detected, this defaults also to the HostExecSpace
 */
using DeviceExecSpace = Kokkos::DefaultExecutionSpace;

/**
 * Corresponding Memory Space to the DeviceExecSpace
 */
using DeviceMemSpace = DeviceExecSpace::memory_space;

/**
 * Default Kokkos Host Execution Space
 */
using HostExecSpace = Kokkos::HostSpace::execution_space;

/**
 * Corresponding Memory Space to the HostExecSpace
 */
using HostMemSpace = HostExecSpace::memory_space;

}  // namespace autopas::utils::kokkos

#endif
