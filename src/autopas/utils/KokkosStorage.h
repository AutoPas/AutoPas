/**
 * @file KokkosStorage.h
 * @author Luis Gall
 * @date 10.12.2025
 */

#pragma once

#include <Kokkos_Core.hpp>

#include "KokkosAoS.h"
#include "KokkosDataLayoutConverter.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas::utils {

  template <class MemSpace, class Particle_T>
  class KokkosStorage {

  public:
    KokkosStorage() {}

    void resize(size_t numParticles) {
      switch (_layout) {
        case DataLayoutOption::aos: {
          storageAoS.resize(numParticles);
          break;
        }
        case DataLayoutOption::soa: {
          storageSoA.resize(numParticles);
          break;
        }
      }
    }

    void addParticle(size_t index, const Particle_T &p) {
      switch (_layout) {
        case DataLayoutOption::aos: {
          storageAoS.getParticle(index) = p;
          break;
        }
        case DataLayoutOption::soa: {
          storageSoA.addParticle(index, p);
          break;
        }
      }
    }

    void convertToSoA(size_t size) {
      constexpr size_t tupleSize = storageSoA.tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      storageSoA.resize(size);
      _converter.convertToSoA(storageAoS, storageSoA, size, I);
    }

    void convertToAoS(size_t size) {
      constexpr size_t tupleSize = storageSoA.tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      storageAoS.resize(size);
      _converter.convertToAoS(storageSoA, storageAoS, size, I);
    }

    template <size_t attribute, bool offset, bool host = false>
    KOKKOS_INLINE_FUNCTION
    constexpr auto& operator() (int i) const {
      switch (_layout) {
        case DataLayoutOption::aos: {
          return storageAoS.template operator()<attribute, offset, host>(i);
        }
        case DataLayoutOption::soa: {
          return storageSoA.template operator()<attribute, offset, host>(i);
        }
        default: {
          // THIS SHOULD NEVER HAPPEN, TODO: log an error
          return storageAoS.template operator()<attribute, offset, host>(i);
        }
      }
    }

    /*
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
    */

    void setLayout(DataLayoutOption newLayout) {
      _layout = newLayout;
    }

    DataLayoutOption getLayout() const {
      return _layout;
    }

    KokkosAoS<MemSpace, Particle_T>& getAoS() {
      return storageAoS;
    }

    const KokkosAoS<MemSpace, Particle_T>& getAoS() const {
      return storageAoS;
    }

    Particle_T::template KokkosSoAArraysType<MemSpace>& getSoA() {
      return storageSoA;
    }

    const Particle_T::template KokkosSoAArraysType<MemSpace>& getSoA() const {
      return storageSoA;
    }

    constexpr static size_t tupleSize() {
      return Particle_T::template KokkosSoAArraysType<MemSpace>::tupleSize();
    }

    size_t size() const {
      switch (_layout) {
        case DataLayoutOption::aos: {
          return storageAoS.size();
        }
        case DataLayoutOption::soa: {
          return storageSoA.size();
        }
        default: {
          return 0;
        }
      }
    }

  private:
    KokkosDataLayoutConverter<Particle_T> _converter {};

    DataLayoutOption _layout {DataLayoutOption::aos};

    KokkosAoS<MemSpace, Particle_T> storageAoS {};
    Particle_T::template KokkosSoAArraysType<MemSpace> storageSoA {};
  };

}