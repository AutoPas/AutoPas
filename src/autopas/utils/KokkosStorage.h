/**
 * @file KokkosStorage.h
 * @author Luis Gall
 * @date 10.12.2025
 */

#pragma once

#include <Kokkos_Core.hpp>

#include "KokkosAoS.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas::utils {

  template <class MemSpace, class Particle_T>
  class KokkosStorage {

  public:
    KokkosStorage() {}

    template <size_t attribute, bool offset>
    KOKKOS_INLINE_FUNCTION
    constexpr auto get(size_t index) {
      switch (_layout) {
        case DataLayoutOption::aos: {
          return storageAoS.template get<attribute, offset>(index);
        }
        case DataLayoutOption::soa: {
          return storageSoA.template get<attribute, offset>(index);
        }
        default: {
          // Should never happen, just to avoid compiler warnings
          return storageAoS.template get<attribute, offset>(index);
        }
      }
    }

    template <size_t attribute, bool offset>
    KOKKOS_INLINE_FUNCTION
    void set(auto value, size_t index) {
      switch (_layout) {
        case DataLayoutOption::aos: {
          storageAoS.template set<attribute, offset>(value, index);
          break;
        }
        case DataLayoutOption::soa: {
          storageSoA.template set<attribute, offset>(value, index);
          break;
        }
      }
    }

    void setLayout(DataLayoutOption newLayout) {
      _layout = newLayout;
    }

  private:
    DataLayoutOption _layout {DataLayoutOption::aos};

    KokkosAoS<MemSpace, Particle_T> storageAoS {};
    Particle_T::template KokkosSoAArraysType<MemSpace> storageSoA {};
  };

}