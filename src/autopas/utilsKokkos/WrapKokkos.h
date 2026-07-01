/**
* @file WrapKokkos.h
 * @author Luis Gall
 * @date 04.11.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#define AUTOPAS_KOKKOS_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#else
#define AUTOPAS_KOKKOS_INLINE_FUNCTION inline
#endif

namespace autopas {

inline int AutoPas_Kokkos_Init(int& argc, char** argv) {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  return 0;
}

inline int AutoPas_Kokkos_Finalize() {
#ifdef AUTOPAS_ENABLE_KOKKOS
  Kokkos::finalize();
#endif
  return 0;
}

}