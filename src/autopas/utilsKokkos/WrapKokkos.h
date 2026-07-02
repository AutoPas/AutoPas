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
 * Wrapper for Kokkos::initialize()
 * @param argc: reference to number of arguments
 * @param argv: argument vector
 */
inline void AutoPas_Kokkos_Init(int &argc, char **argv) {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
}

/**
 * Wrapper for Kokkos::finalize()
 */
inline void AutoPas_Kokkos_Finalize() {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::finalize();
#endif
}

}  // namespace autopas