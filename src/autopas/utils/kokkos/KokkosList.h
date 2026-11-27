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
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::utils::kokkos {

/**
 * Kokkos List class that wraps different data layouts into one class and delegates requests to the currently active
 * format
 * @tparam Particle_T The particle type to store
 */
template <typename Particle_T>
class KokkosList {
 public:
  KokkosList() {}

  /**
   * Constructor for KokkosList for a specific data layout and capacity
   * @param activeDataLayout Data Layout for the initial list
   * @param size initial capacity of the list
   */
  KokkosList(DataLayoutOption activeDataLayout, size_t size) : _activeDataLayout(activeDataLayout) { resize(size); }

  /**
   * Operator to access the data of the requested particle attribute at the given position.
   * Delegates the operation to the currently active data layout
   * @tparam Attribute Attribute name
   * @tparam useDeviceView Flag to choose whether Kokkos should access the device-side memory
   * @param index particle index in the list
   * @return value of the particle attribute in the list at position index
   */
  template <typename Particle_T::AttributeNames Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto &operator()(size_t index) const {
    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        return _aos.template operator()<Attribute, useDeviceView>(index);
      }
      case DataLayoutOption::Value::soa: {
        return _soa.template operator()<Attribute, useDeviceView>(index);
      }
    }
    Kokkos::abort("KokkosList::operator(): invalid data layout");
    // Kokkos::abort is not [[noreturn]] on NVHPC, SYCL device passes, and OpenACC; unreachable otherwise.
    return _soa.template operator()<Attribute, useDeviceView>(index);
  }

  /**
   * Sets the capacity of the list to the given size. If size < the current size, this throws an exception.
   * Delegates the operation to the currently active data layout.
   * @param size new capacity of the list
   */
  void resize(size_t size) {
    if (size < _size) {
      utils::ExceptionHandler::exception(
          "KokkosList::resize(): invalid size {}. The requested size is smaller than the current size. This would lead "
          "to particles being deleted silently, which is not valid.",
          size);
    }

    _capacity = size;
    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        _aos.resize(size);
        return;
      }
      case DataLayoutOption::Value::soa: {
        _soa.resize(size);
        return;
      }
    }
    utils::ExceptionHandler::exception("KokkosList::resize(): invalid data layout {}",
                                       static_cast<int>(_activeDataLayout));
  }

  /**
   * Marks one attribute as modified in the given space. Delegates the operation to the currently active data layout.
   * The AoS can only be marked as a whole, hence the attribute is ignored for AoS.
   * @tparam Space
   * @tparam Attribute
   */
  template <typename Space, typename Particle_T::AttributeNames Attribute>
  void modify() {
    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        _aos.template modify<Space>();
        return;
      }
      case DataLayoutOption::Value::soa: {
        _soa.template modify<Space, Attribute>();
        return;
      }
    }
    utils::ExceptionHandler::exception("KokkosList::modify(): invalid data layout {}",
                                       static_cast<int>(_activeDataLayout));
  }

  /**
   * Syncs one attribute to the given space. Delegates the operation to the currently active data layout.
   * The AoS can only be synced as a whole, hence the attribute is ignored for AoS.
   * @tparam Space
   * @tparam Attribute
   */
  template <typename Space, typename Particle_T::AttributeNames Attribute>
  void sync() {
    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        _aos.template sync<Space>();
        return;
      }
      case DataLayoutOption::Value::soa: {
        _soa.template sync<Space, Attribute>();
        return;
      }
    }
    utils::ExceptionHandler::exception("KokkosList::sync(): invalid data layout {}",
                                       static_cast<int>(_activeDataLayout));
  }

  /**
   * Syncs all attributes to the given space.
   * @tparam Space
   */
  template <typename Space>
  void syncAll() {
    _aos.template sync<Space>();
    _soa.template syncAll<Space>();
  }

  /**
   * Marks all attributes as modified in the given space.
   * @tparam Space
   */
  template <typename Space>
  void modifyAll() {
    _aos.template modify<Space>();
    _soa.template modifyAll<Space>();
  }

  /**
   * Getter for the number of particles
   * @return number of particles in the list
   */
  size_t size() const { return _size; }

  /**
   * Getter for the capacity of the list
   * @return capacity of the list
   */
  size_t getCapacity() const { return _capacity; }

  /**
   * Adds the given particle to the list.
   * If the capacity is exceeded, a resize to double the capacity of the list is executed first.
   * Since particles can only be added on the host, this operation syncs the particle data to the Host and marks it as
   * modified on the Host. Delegates the operation to the currently active data layout
   * @param p Particle object to add to the list
   */
  void addParticle(const Particle_T &p) {
    if (_size == _capacity) {
      size_t newCapacity = _capacity == 0 ? 2 : _capacity * 2;
      resize(newCapacity);
      _capacity = newCapacity;
    }

    syncAll<HostExecSpace>();
    modifyAll<HostExecSpace>();

    switch (_activeDataLayout) {
      case DataLayoutOption::Value::aos: {
        return _aos.addParticle(_size++, p);
      }
      case DataLayoutOption::Value::soa: {
        return _soa.addParticle(_size++, p);
      }
    }
    utils::ExceptionHandler::exception("KokkosList::addParticle(): invalid data layout {}",
                                       static_cast<int>(_activeDataLayout));
  }

  /**
   * Getter for the currently active data layout
   * @return currently active data layout
   */
  DataLayoutOption getActiveDataLayout() const { return _activeDataLayout; }

  /**
   * Getter for the KokkosAoS list
   * @return KokkosAoS list
   */
  const auto &getAoS() const { return _aos; }

  /**
   *
   * @return KokkosSoA list
   */
  const auto &getSoA() const { return _soa; }

 private:
  /**
   * Number of particles in the list
   */
  size_t _size{0};

  /**
   * Capacity of the list
   */
  size_t _capacity{0};

  /**
   * Currently active data layout
   */
  DataLayoutOption _activeDataLayout{DataLayoutOption::soa};

  /**
   * AoS list
   */
  KokkosAoS<Particle_T> _aos{};

  /**
   * SoA list
   */
  Particle_T::KokkosSoAArraysType _soa{};
};

}  // namespace autopas::utils::kokkos

#endif