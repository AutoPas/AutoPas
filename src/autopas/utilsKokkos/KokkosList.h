/**
 * @file KokkosList.h
 * @author Luis Gall
 * @date 10.12.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>

namespace autopas::utilsKokkos {

template <typename Particle_T>
class KokkosList {

  template <size_t Attribute>
  KOKKOS_INLINE_FUNCTION auto& operator() (size_t index) const {

  }

  void shrinkToFit() {
    // todo: adapt capacity to size
  }

  void reserve(size_t size) {

  }

  size_t size() const {
    return _size;
  }

  void addParticle(Particle_T p) {

  }

private:
  size_t _size {0};
  size_t _capacity {0};

};

} // utilsKokkos

#endif