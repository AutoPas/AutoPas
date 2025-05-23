/**
 * Contains wrapper for the Kokkos functions
 * @file WrapKokkos.h
 * @author schuhmaj
 * @date 23 May 2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS
#include "Kokkos_Core.hpp"
#endif

namespace autopas::kokkos {

/**
 * Simple wrapper to Kokkos initialize.
 * This function needs to be called before anything else.
 * If MPI is enabled. This function needs to be called after MPI_Init!
 * @param argc the argument count (modifiable as Kokkos might remove arguments)
 * @param argv the arguments (modifiable as Kokkos might remove arguments)
 */
inline void initialize(int &argc, char *argv[]) {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
}

/**
 * Simple wrapper for Kokkos finalize.
 * This function needs to be called after everything else.
 * If MPI is enabled, this function needs to be called before MPI_Finalize!
 */
inline void finalize() {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::finalize();
#endif
}

}  // namespace autopas::kokkos