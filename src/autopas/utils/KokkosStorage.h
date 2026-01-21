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

  template <class Particle_T>
  class KokkosStorage {

  public:
    KokkosStorage() {}

    KokkosStorage(const KokkosStorage& other) {
      _layout = other.getLayout();

      switch (_layout) {
        case DataLayoutOption::aos: {
          storageAoS = other.getAoS();
          break;
        }
        case DataLayoutOption::soa: {
          storageSoA = other.getSoA();
          break;
        }
      }
    }

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
          storageAoS.addParticle(index, p);
          break;
        }
        case DataLayoutOption::soa: {
          sync<Kokkos::HostSpace::execution_space>();
          storageSoA.addParticle(index, p);
          markModified<Kokkos::HostSpace::execution_space>();
          break;
        }
      }
    }

    void convertToSoA(size_t size) {
      constexpr size_t tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      storageSoA.resize(size);
      _converter.convertToSoA(storageAoS, storageSoA, size, I);
    }

    void convertToAoS(size_t size) {
      constexpr size_t tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      storageAoS.resize(size);
      _converter.convertToAoS(storageSoA, storageAoS, size, I);
    }

    template <size_t attribute, bool offset, bool host = false>
    KOKKOS_INLINE_FUNCTION
    constexpr auto& operator() (int i) const {
      switch (_layout) {
        case DataLayoutOption::aos: {
          return storageAoS.template operator()<attribute, offset>(i);
        }
        case DataLayoutOption::soa: {
          return storageSoA.template operator()<attribute, offset, host>(i);
        }
        default: {
          // THIS SHOULD NEVER HAPPEN, TODO: log an error
          return storageAoS.template operator()<attribute, offset>(i);
        }
      }
    }

    Particle_T& getParticle(int i) const {
      return storageAoS.getParticle(i);
    }

    template <typename Target>
    void sync() {
      constexpr auto tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      storageSoA.template syncAll<Target>(I);
    }

    template <typename Target>
    void markModified() {
      constexpr auto tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      storageSoA.template markAllModified<Target>(I);
    }

    void setLayout(DataLayoutOption newLayout) {
      _layout = newLayout;
    }

    DataLayoutOption getLayout() const {
      return _layout;
    }

    KokkosAoS<Particle_T>& getAoS() {
      return storageAoS;
    }

    const KokkosAoS<Particle_T>& getAoS() const {
      return storageAoS;
    }

    Particle_T::KokkosSoAArraysType& getSoA() {
      return storageSoA;
    }

    const Particle_T::KokkosSoAArraysType& getSoA() const {
      return storageSoA;
    }

    constexpr static size_t tupleSize() {
      return Particle_T::KokkosSoAArraysType::tupleSize();
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

    KokkosAoS<Particle_T> storageAoS {};
    Particle_T::KokkosSoAArraysType storageSoA {};
  };

}