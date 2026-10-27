/**
 * @file ReferenceCellIterator.h
 * @author S. Newcome
 * @date 15.10.2026
 */

#pragma once

#include <iterator>
#include <type_traits>
#include <vector>

namespace autopas {

/**
 * Iterator for "references"-based cells, i.e. where the cell itself does not contain the particles, but all particles
 * are instead stored in a global std::vector, and a cell's particles occupy a contiguous segment of that vector. It
 * provides the same interface as CellIterator.
 *
 * @note The iterator stores the address of the particle vector's buffer, so adding particles to that vector while
 * iterating invalidates the iterator (this typically holds for all containers).
 *
 * @tparam Particle_T Type of the particles.
 * @tparam modifiable If false, this is a const iterator, meaning, the underlying particle can not be modified.
 */
template <class Particle_T, bool modifiable>
class ReferenceCellIterator {
 public:
  /**
   * Type of the wrapped iterator over the indices.
   */
  using IteratorType =
      std::conditional_t<modifiable, std::vector<size_t>::iterator, std::vector<size_t>::const_iterator>;

  /**
   * Type trait for std::iterator_traits, to enable use of standard library algorithms. Must be of the name expected by
   * the standard library, hence not camel case like the rest of the project.
   */
  using difference_type = IteratorType::difference_type;
  /**
   * Type trait for std::iterator_traits, to enable use of standard library algorithms. The particle type rather than
   * the stored index type, and non-const even for const iterators. Must be of the name expected by the standard
   * library, hence not camel case like the rest of the project.
   */
  using value_type = Particle_T;
  /**
   * Type trait for std::iterator_traits, to enable use of standard library algorithms.
   */
  using pointer = std::conditional_t<modifiable, Particle_T *, const Particle_T *>;
  /**
   * Type trait for std::iterator_traits, to enable use of standard library algorithms.
   */
  using reference = std::conditional_t<modifiable, Particle_T &, const Particle_T &>;
  /**
   * Type trait for std::iterator_traits, to enable use of standard library algorithms. Must be of the name expected by
   * the standard library, hence not camel case like the rest of the project.
   */
  using iterator_category = std::bidirectional_iterator_tag;

  /**
   * Constructor
   * @param particles Start of the buffer of the vector the indices refer to.
   * @param iterator Iterator over the indices.
   */
  ReferenceCellIterator(pointer particles, IteratorType iterator) : _particles(particles), _iterator(iterator) {}

  /**
   * Dereference operator.
   * @return Reference to the current particle.
   */
  reference operator*() const { return _particles[*_iterator]; }

  /**
   * Dereference operator.
   * @return Pointer to the current particle.
   */
  pointer operator->() const { return &operator*(); }

  /**
   * Increment the iterator.
   * @return *this
   */
  ReferenceCellIterator &operator++() {
    ++_iterator;
    return *this;
  }

  /**
   * Decrement the iterator.
   * @return *this
   */
  ReferenceCellIterator &operator--() {
    --_iterator;
    return *this;
  }

  /**
   * Equality operator
   * @param rhs
   * @return
   */
  bool operator==(const ReferenceCellIterator &rhs) const { return _iterator == rhs._iterator; }

  /**
   * Not equality operator
   * @param rhs
   * @return
   */
  bool operator!=(const ReferenceCellIterator &rhs) const { return not(*this == rhs); }

  /**
   * Distance between two iterators.
   * @param rhs
   * @return Number of elements between two iterators.
   */
  difference_type operator-(const ReferenceCellIterator &rhs) const { return _iterator - rhs._iterator; }

  /**
   * Comparison operator.
   * @param rhs
   * @return True if this is before (in -- direction) rhs.
   */
  bool operator<(const ReferenceCellIterator &rhs) const { return _iterator < rhs._iterator; }

 private:
  /**
   * Start of the buffer of the vector the indices refer to.
   */
  pointer _particles;

  /**
   * Actual underlying iterator over the indices.
   */
  IteratorType _iterator;
};

}  // namespace autopas