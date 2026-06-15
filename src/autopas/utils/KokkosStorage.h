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
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/inBox.h"

namespace autopas::utils {

  template <class Particle_T>
  class KokkosStorage {

  public:
    // TODO: think about deleting this constructor
    KokkosStorage() {}

    KokkosStorage(DataLayoutOption layout, size_t numParticles) : _layout(layout) {
      resize(numParticles);
    }

    KokkosStorage(const KokkosStorage& other) {
      _layout = other.getLayout();
      _currentSize = other.size();
      _capacity = other.getCapacity();

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

    void realloc(size_t numParticles) {
      switch (_layout) {
        case DataLayoutOption::aos: {
          storageAoS.realloc(numParticles);
          break;
        }
        case DataLayoutOption::soa: {
          storageSoA.realloc(numParticles);
          break;
        }
      }

      _capacity = numParticles;
      _currentSize = 0; // because realloc deletes all content
    }

    void resize(size_t numParticles) {
      // TODO: check if numParticles < currentSize and decide how to handle that

      switch (_layout) {
        case DataLayoutOption::aos: {
          storageAoS.resize(numParticles);
          _aosDirty = true;
          break;
        }
        case DataLayoutOption::soa: {
          storageSoA.resize(numParticles);
          _soaDirty = true;
          break;
        }
      }

      _currentSize = std::min(_currentSize, numParticles);
      _capacity = numParticles;
    }

    void clear () {
      storageAoS.realloc(0);
      storageSoA.realloc(0);

      _capacity = 0;
      _currentSize = 0;

      _soaDirty = false;
      _aosDirty = false;
    }

    // TODO: mark somewhere that this is not thread safe (!)
    void addParticle(const Particle_T &p) {

      if (_currentSize == _capacity) {
        resize(_capacity+10);
      }

      size_t index = _currentSize++;

      switch (_layout) {
        case DataLayoutOption::aos: {
          syncSoAToAoS();
          storageAoS.addParticle(index, p);
          _aosDirty = true;
          break;
        }
        case DataLayoutOption::soa: {
          syncAoSToSoA();
          sync<Kokkos::HostSpace::execution_space>();
          storageSoA.addParticle(index, p);
          markModified<Kokkos::HostSpace::execution_space>();
          _soaDirty = true;
          break;
        }
      }
    }

    void syncAoSToSoA() {
      if (_aosDirty) {
        constexpr size_t tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
        constexpr auto I = std::make_index_sequence<tupleSize>();

        const auto size = storageAoS.size();
        storageSoA.resize(size);
        _converter.convertToSoA(storageAoS, storageSoA, size, I);
        _aosDirty = false;
      }
    }

    void syncSoAToAoS() {
      if (_soaDirty) {
        constexpr size_t tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
        constexpr auto I = std::make_index_sequence<tupleSize>();

        const auto size = storageSoA.size();
        storageAoS.resize(size);
        _converter.convertToAoS(storageSoA, storageAoS, size, I);
        _soaDirty = false;
      }
    }

    template <size_t attribute, bool offset, bool host = false>
    KOKKOS_INLINE_FUNCTION
    constexpr auto& operator() (int i) const {
      switch (_layout) {
        case DataLayoutOption::aos: {
          return storageAoS.template operator()<attribute>(i);
        }
        case DataLayoutOption::soa: {
          return storageSoA.template operator()<attribute, offset, host>(i);
        }
        default: {
          // THIS SHOULD NEVER HAPPEN, TODO: log an error
          return storageAoS.template operator()<attribute>(i);
        }
      }
    }

    template <bool host>
    void copyParticle(int targetIndex, const KokkosStorage<Particle_T>& otherStorage, int sourceIndex) const {
      constexpr auto tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      this->template copyParticleImpl<host>(targetIndex, otherStorage, sourceIndex, I);
    }

    template <bool regionIter, bool host, typename T>
    KOKKOS_INLINE_FUNCTION
    bool fulfillsIteratorRequirements(int index, autopas::options::IteratorBehavior behavior, const Kokkos::Array<T, 3>& lowerCorner, const Kokkos::Array<T, 3>& upperCorner) const {

      if constexpr (regionIter) {
        std::array<T, 3> positions {
          operator()<Particle_T::AttributeNames::posX, true, host>(index),
          operator()<Particle_T::AttributeNames::posY, true, host>(index),
          operator()<Particle_T::AttributeNames::posZ, true, host>(index),
        };

        // TODO: write own version of inBox maybe in autopas::kokkosUtils namespace

        /*
        if (not autopas::utils::inBox(positions, lowerCorner, upperCorner)) {
          return false;
        }
        */
      }

      auto ownershipState = operator()<Particle_T::AttributeNames::ownershipState, true, host>(index);

      // TODO: this will require checks for edge cases and sync with the changes for dummy particles
      return static_cast<unsigned int>(ownershipState) & static_cast<unsigned int>(behavior);
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
      // TODO: when AoS is also allowed on the GPU, we need a switch here too
      constexpr auto tupleSize = Particle_T::KokkosSoAArraysType::tupleSize();
      constexpr auto I = std::make_index_sequence<tupleSize>();

      storageSoA.template markAllModified<Target>(I);
    }

    // TODO: rethink this interface and when this should be used
    void markLayoutModified(DataLayoutOption layout) {
      if (layout == DataLayoutOption::soa) {
        _soaDirty = true;
      } else if (layout == DataLayoutOption::aos) {
        _aosDirty = true;
      }
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

    // TODO: make it very (!) clear that this function is only to be used in exceptional cases
    void overrideSize(size_t newSize) {
      _currentSize = newSize;
    }

    void overrideCapacity(size_t newCapacity) {
      _capacity = newCapacity;
    }

    size_t size() const {
      return _currentSize;
    }

    size_t getCapacity() const {
      return _capacity;
    }

  private:

    template <bool host, std::size_t... I>
    void copyParticleImpl (int targetIndex, const KokkosStorage<Particle_T>& otherStorage, int sourceIndex, std::index_sequence<I...>) const {

      switch (_layout) {
        case DataLayoutOption::aos: {
          ((this->storageAoS.template operator()<I+1>(targetIndex) = otherStorage.getAoS().template operator()<I+1>(sourceIndex)), ...);
          break;
        }
        case DataLayoutOption::soa: {
          ((this->storageSoA.template operator()<I, false, host>(targetIndex) = otherStorage.getSoA().template operator()<I, false, host>(sourceIndex)), ...);
          break;
        }
      }
    }

    KokkosDataLayoutConverter<Particle_T> _converter {};

    DataLayoutOption _layout {DataLayoutOption::aos};

    KokkosAoS<Particle_T> storageAoS {};
    Particle_T::KokkosSoAArraysType storageSoA {};

    // TODO: maybe think of a better concept for this in the future
    bool _soaDirty {false};
    bool _aosDirty {false};

    size_t _capacity {0};
    size_t _currentSize {0};
  };

}