/**
 * @file WrapKokkos.h
 * @author Luis Gall
 * @date 04.11.2025
 */

#pragma once

#include <Kokkos_Core.hpp>

namespace autopas {

inline int AutoPas_Kokkos_Init(int& argc, char** argv) {
  Kokkos::initialize(argc, argv);
  return 0;
}

inline int AutoPas_Kokkos_Finalize() {
  Kokkos::finalize();
  return 0;
}



}