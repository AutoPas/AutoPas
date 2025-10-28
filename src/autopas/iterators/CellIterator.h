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
 * @note This class is necessary because some cells store `Particle_T` and others `Particle_T *`.
 * @tparam StorageType Type of the cell, used to derive the particle type.
 * @tparam modifiable If false, this is a const iterator, meaning, the underlying particle can not be modified.
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
   * Type trait to be compatible with std::iterator.
   */
  using difference_type = typename IteratorType::difference_type;
  /**
   * Type trait to be compatible with std::iterator.
   */
  using value_type = typename IteratorType::value_type;
  /**
   * Type trait to be compatible with std::iterator.
   */
  using pointer = typename IteratorType::pointer;
  /**
   * Type trait to be compatible with std::iterator.
   */
  using reference = typename IteratorType::reference;
  /**
   * Type trait to be compatible with std::iterator.
   */
  using iterator_category = typename IteratorType::iterator_category;

  /**
   * Constructor
   * @param iterator
   */
  explicit CellIterator(IteratorType iterator) : iterator(iterator) {}

  /**
   * Dereference operator.
   * @return Reference to the current particle.
   */
  inline std::conditional_t<modifiable, ParticleType &, const ParticleType &> operator*() const {
    using RetType = ParticleType &;
    // depending on the cell the iterator points to a particle or pointer
    // resolve this at compile time
    if constexpr (std::is_same_v<std::remove_const_t<std::remove_reference_t<RetType>>,
                                 std::remove_const_t<std::remove_reference_t<decltype(*iterator)>>>) {
      return *iterator;
    } else if constexpr (std::is_same_v<std::remove_const_t<std::remove_reference_t<RetType>>,
                                        std::remove_const_t<std::remove_reference_t<decltype(**iterator)>>>) {
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
  inline std::conditional_t<modifiable, ParticleType *, const ParticleType *> operator->() const {
    return &operator*();
  }

  /**
   * Increment the iterator.
   * @return *this
   */
  inline CellIterator<StorageType, modifiable> &operator++() {
    ++iterator;
    return *this;
  }

  /**
   * Decrement the iterator.
   * @return *this
   */
  inline CellIterator<StorageType, modifiable> &operator--() {
    --iterator;
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

  /**
   * Distance between two iterators.
   * @param rhs
   * @return Number of elements between two iterators.
   */
  difference_type operator-(const CellIterator &rhs) const { return iterator - rhs.iterator; }

  /**
   * Comparison operator.
   * @param rhs
   * @return True if this is before (in -- direction) rhs.
   */
  bool operator<(const CellIterator &rhs) const { return iterator < rhs; }

 private:
  /**
   * Actual underlying iterator.
   */
  IteratorType iterator;
};
