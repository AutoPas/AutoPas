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

  template <size_t Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto& operator() (size_t index) const {
    if constexpr(useDeviceView) {
      return particles.view_device()(index).template operator()<static_cast<Particle_T::AttributeNames>(Attribute)>();
    }
    return particles.view_host()(index).template operator()<static_cast<Particle_T::AttributeNames>(Attribute)>();
  }

  void resize(size_t size) {
    particles.resize(Kokkos::WithoutInitializing, size);
  }

  template <bool useDeviceView = true>
  void addParticle(size_t index, Particle_T& particle) {
    if constexpr(useDeviceView) {
      particles.view_device()(index) = particle;
    } else {
      particles.view_host()(index) = particle;
    }
  }

private:
  Kokkos::DualView<Particle_T, DeviceExecSpace::device_type> particles {};
};

}

#endif