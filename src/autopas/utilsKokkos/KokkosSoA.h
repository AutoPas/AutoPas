/**
* @file KokkosSoA.h
 * @author Luis Gall
 * @date 13.11.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>
#include "Kokkos_DualView.hpp"

#include "autopas/utilsKokkos/KokkosSpaces.h"

namespace autopas::utilsKokkos {

template <typename ... AttributeTypes>
class KokkosSoA {


private:
  Kokkos::DualView<AttributeTypes..., DeviceExecSpace::device_type> particleData {};
};

} // utilsKokkos

#endif