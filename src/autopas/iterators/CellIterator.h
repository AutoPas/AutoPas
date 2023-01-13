/**
 * @file CellIterator.h
 * @author F. Gratl
 * @date 12.01.23
 */

#pragma once

#include <type_traits>
#include <vector>

/**
 * Wraps the iterator of arbitrary cells to provide a common interface.
 * @note This class is necessary because some cells store `Particle` and others `Particle *`.
 * @tparam CellType Type of the cell, used to derive the particle and storage type.
 * @tparam modifiable
 */
template <class StorageType, bool modifiable>
class CellIterator {
 public:
  /**
   * Type of the wrapped iterator.
   */
  using IteratorType =
      std::conditional_t<modifiable, typename StorageType::iterator, typename StorageType::const_iterator>;

  /**
   * Type of the Particles in the storage. Should always be the actual type and not a pointer.
   */
  using ParticleType = std::remove_pointer_t<typename StorageType::value_type>;

  /**
   * Constructor
   * @param iterator
   */
  explicit CellIterator(IteratorType iterator) : iterator(iterator) {}

  /**
   * Dereference operator.
   * @return Reference to the current particle.
   */
  inline ParticleType &operator*() const {
    using RetType = ParticleType &;
    // depending on the cell the iterator points to a particle or pointer
    // resolve this at compile time
    if constexpr (std::is_same_v<RetType, decltype(*iterator)>) {
      return *iterator;
    } else if constexpr (std::is_same_v<RetType, decltype(**iterator)>) {
      return **iterator;
    } else {
      // trigger compile error that prints the type of *iterator
      return decltype(*iterator)::thisMemberDoesntExist;
    }
  }

  /**
   * Dereference operator.
   * @return Pointer to the current particle.
   */
  inline ParticleType *operator->() const { return &operator*(); }

  /**
   * Increments the iterator.
   * @return *this
   */
  inline CellIterator<StorageType, modifiable> &operator++() {
    ++iterator;
    return *this;
  }

  /**
   * Equality operator
   * @param rhs
   * @return
   */
  bool operator==(const CellIterator &rhs) const { return iterator == rhs.iterator; }
  /**
   * Not equality operator
   * @param rhs
   * @return
   */
  bool operator!=(const CellIterator &rhs) const { return not(*this == rhs); }

 private:
  IteratorType iterator;
};