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
 public:
  KokkosList() : _activeDataLayout(DataLayoutOption::soa) {}

  KokkosList(DataLayoutOption activeDataLayout, size_t size) : _activeDataLayout(activeDataLayout) { resize(size); }

  template <size_t Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto &operator()(size_t index) const {
    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        return _aos.template operator()<Attribute, useDeviceView>(index);
      }
      case DataLayoutOption::Value::soa: {
        return _soa.template operator()<Attribute-1, useDeviceView>(index);
      }
      default:
        break;  // todo: throw exception although this should currently not happen
    }
  }

  void resize(size_t size) {
    if (size < _size) {
      // todo: warn about that the new size is smaller than the current size and that particle data is lost by this
      // action todo: maybe even throw and exception since this can lead to bugs/undesired behavior when people are not
      // aware of that
      _size = size;
    }

    _capacity = size;
    _soa.resize(size);
    _aos.resize(size);
  }

  template <typename Space, std::size_t I>
  void modify() {
    _soa.template modify<Space>(I);
    _aos.template modify<Space>();
  }

  template <typename Space, std::size_t I>
  void sync() {
    _soa.template sync<Space>(I);
    _aos.template sync<Space>();
  }

  template <typename Space>
  void syncAll() {
    _aos.template sync<Space>();
    _soa.template syncAll<Space>();
  }

  template <typename Space>
  void modifyAll() {
    _aos.template modify<Space>();
    _soa.template modifyAll<Space>();
  }

  size_t size() const { return _size; }

  size_t getCapacity() const { return _capacity; }

  void addParticle(const Particle_T &p) {

    syncAll<HostExecSpace>();
    modifyAll<HostExecSpace>();

    if (_size == _capacity) {
      size_t newCapacity = _capacity == 0 ? 2 : _capacity * 2;
      resize(newCapacity);
      _capacity = newCapacity;
    }

    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        return _aos.addParticle(_size++, p);
      }
      case DataLayoutOption::Value::soa: {
        return _soa.addParticle(_size++, p);
      }
      default:
        break;  // todo: throw exception although this should currently not happen
    }
  }

  DataLayoutOption getActiveDataLayout() const { return _activeDataLayout; }

  const auto &getAoS() const { return _aos; }

  const auto &getSoA() const { return _soa; }

 private:

  size_t _size{0};
  size_t _capacity{0};

  DataLayoutOption _activeDataLayout{};

  KokkosAoS<Particle_T> _aos{};
  Particle_T::KokkosSoAArraysType _soa{};
};

}  // namespace autopas::utilsKokkos

#endif