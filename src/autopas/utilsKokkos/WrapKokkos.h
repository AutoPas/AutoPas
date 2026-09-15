/**
 * @file WrapKokkos.h
 * @author Luis Gall
 * @date 04.11.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

namespace autopas {

/**
 * Wrapper for Kokkos::initialize().
 *
 * Should be called once at the beginning of the program, before any Kokkos-based container or
 * traversal is used. If AutoPas was compiled without Kokkos support, this is a no-op.
 *
 * @param argc: reference to number of arguments
 * @param argv: argument vector
 */
inline void AutoPas_Kokkos_Init([[maybe_unused]] int &argc, [[maybe_unused]] char **argv) {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
}

/**
 * Wrapper for Kokkos::finalize().
 *
 * Should be called exactly once at the end of the program, after all Kokkos work has finished.
 * If AutoPas was compiled without Kokkos support, this is a no-op.
 */
inline void AutoPas_Kokkos_Finalize() {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::finalize();
#endif
}

}  // namespace autopas
