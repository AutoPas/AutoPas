/**
 * @file KokkosList.h
 * @author Luis Gall
 * @date 10.12.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>

#include "KokkosAoS.h"
#include "KokkosSoA.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas::utilsKokkos {

template <typename Particle_T>
class KokkosList {

  template <size_t Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto& operator() (size_t index) const {
    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        return _soa.template operator()<Attribute, useDeviceView>(index);
      }
      case DataLayoutOption::Value::soa: {
        return _soa.template operator()<Attribute, useDeviceView>(index);
      }
      default:
        break; // todo: throw exception although this should currently not happen
    }
  }

  void resize(size_t size) {
    if (size < _size) {
      // todo: warn about that the new size is smaller than the current size and that particle data is lost by this action
      // todo: maybe even throw and exception since this can lead to bugs/undesired behavior when people are not aware of that
      _size = size;
    }

    _capacity = size;
    _soa.resize(size);
    _aos.resize(size);
  }

  size_t size() const {
    return _size;
  }

  void addParticle(Particle_T& p) {
    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        return _soa.addParticle<true>(p);
      }
      case DataLayoutOption::Value::soa: {
        return _soa.addParticle<true>(p);
      }
      default:
        break; // todo: throw exception although this should currently not happen
    }
  }

private:
  size_t _size {0};
  size_t _capacity {0};

  DataLayoutOption _activeDataLayout {};

  KokkosAoS<Particle_T> _aos {};
  KokkosSoA<typename Particle_T::KokkosSoAArraysType> _soa {};

};

} // utilsKokkos

#endif