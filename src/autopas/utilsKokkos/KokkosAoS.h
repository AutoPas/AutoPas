/**
 * @file KokkosAoS.h
 * @author Luis Gall
 * @date 27.11.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>
#include <type_traits>

#include "Kokkos_DualView.hpp"
#include "autopas/utilsKokkos/KokkosSpaces.h"

namespace autopas::utilsKokkos {

template <typename Particle_T>
class KokkosAoS {
  static_assert(std::is_trivially_copyable_v<Particle_T>,
                "KokkosAoS requires trivially copyable particles: elements are allocated WithoutInitializing and "
                "deep-copied between memory spaces, so e.g. virtual functions are not allowed.");

 public:
  template <typename Particle_T::AttributeNames Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto &operator()(size_t index) const {
    if constexpr (useDeviceView) {
      return _particles.view_device()(index).template operator()<Attribute>();
    } else {
      return _particles.view_host()(index).template operator()<Attribute>();
    }
  }

  void resize(size_t size) { _particles.resize(Kokkos::WithoutInitializing, size); }

  void addParticle(size_t index, const Particle_T &particle) { _particles.view_host()(index) = particle; }

  const auto &getView() const { return _particles; }

  template <typename Space>
  void modify() {
    _particles.template modify<Space>();
  }

  template <typename Space>
  void sync() {
    _particles.template sync<Space>();
  }

 private:
  Kokkos::DualView<Particle_T *, DeviceExecSpace::device_type> _particles{};
};

}  // namespace autopas::utilsKokkos

#endif