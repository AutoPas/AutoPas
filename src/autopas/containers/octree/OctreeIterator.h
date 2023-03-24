/**
 * @file OctreeIterator.h
 * @author F. Gratl
 * @date 24.03.23
 */

#pragma once

#include "OctreeNodeWrapper.h"
#include "autopas/iterators/CellIterator.h"

namespace autopas {

/**
 * Wrapper for the cell iterator for Octrees.
 *
 * @tparam StorageType
 * @tparam modifiable
 */
template <class StorageType, bool modifiable>
class OctreeIterator {
 private:
  /**
   * id of the currently iterated cell in their parent's children array
   */
  int _childId;

  /**
   * Iterator for the current cell.
   */
  CellIterator<StorageType, modifiable> _cellIterator;

 public:
  using ParticleType = typename decltype(_cellIterator)::ParticleType;
  using IteratorType = typename decltype(_cellIterator)::IteratorType;

 private:
  /**
   * End of the current cell.
   * @note needs to be defined after `using ParticleType` which needs to be defined after `_cellIterator`.
   */
  OctreeNodeWrapper<ParticleType> *_currentNode;

 public:
  /**
   * Dereference operator.
   * @return Reference to the current particle.
   */
  inline std::conditional_t<modifiable, ParticleType &, const ParticleType &> operator*() const {
    return *_cellIterator;
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
  inline OctreeIterator<StorageType, modifiable> &operator++() {
    ++_cellIterator;
    if (_cellIterator == _currentNode->end()) {
      const auto parent = _currentNode->getRaw()->getParent();
      if (_childId >= 7) {
        // TODO TRAVERSE UP
      }
    }
    return *this;
  }

  /**
   * Equality operator
   * @param rhs
   * @return
   */
  bool operator==(const OctreeIterator &rhs) const { return _cellIterator == rhs._cellIterator; }
  /**
   * Not equality operator
   * @param rhs
   * @return
   */
  bool operator!=(const OctreeIterator &rhs) const { return not(*this == rhs); }
};
}  // namespace autopas