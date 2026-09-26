/**
* @file KokkosAoS.h
 * @author Luis Gall
 * @date 27.11.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>
#include "Kokkos_DualView.hpp"

#include "autopas/utilsKokkos/KokkosSpaces.h"

namespace autopas::utilsKokkos {

template <typename Particle_T>
class KokkosAoS {


private:
  Kokkos::DualView<Particle_T, DeviceExecSpace::device_type> particles {};
};

}

#endif