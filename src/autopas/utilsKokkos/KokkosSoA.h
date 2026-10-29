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
public:
  using Type = std::tuple<Kokkos::DualView<AttributeTypes, DeviceExecSpace::device_type>...>;

  template <size_t Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto& operator() (size_t index) const {
    if constexpr(useDeviceView) {
      return std::get<Attribute>(particleData).view_device()(index);
    }
    return std::get<Attribute>(particleData).view_host()(index);
  }

  void resize(size_t size) {
    constexpr auto tupleSize = std::tuple_size<decltype(particleData)>::value;
    constexpr auto I = std::make_index_sequence<tupleSize>{};

    resizeImpl(size, I);
  }

  template <class Particle_T, bool useDeviceView>
  void addParticle(size_t index, Particle_T& particle) {
    constexpr auto tupleSize = std::tuple_size<decltype(particleData)>::value;
    constexpr auto I = std::make_index_sequence<tupleSize>{};

    addParticleImpl<Particle_T, useDeviceView>(index, particle, I);
  }

private:
  template <std::size_t... I>
  void resizeImpl(size_t size, std::index_sequence<I...>) {
    (std::get<I>(particleData).resize(Kokkos::WithoutInitializing, size), ...);
  }

  template <class Particle_T, bool useDeviceView, std::size_t... I>
  void addParticleImpl(size_t index, Particle_T& particle, std::index_sequence<I...>) {
    // I+1 is required since we do not want to store the Particle Pointer which is the result of particle.get<0>()
    ((operator()<I, useDeviceView>(index) = particle.template get<static_cast<Particle_T::AttributeNames>(I+1)>()), ...);
  }

  Type particleData {};
};

} // utilsKokkos

#endif