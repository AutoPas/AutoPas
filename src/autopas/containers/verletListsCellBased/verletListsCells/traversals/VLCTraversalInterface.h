/**
 * @file VLCTraversalInterface.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

namespace autopas {

/**
 * This class provides the Traversal Interface for the verlet lists cells container.
 *
 * This class handles traversals through the cell structures with neighbor lists.
 * Derived classes handle the order through which the cells are traversed.
 *
 * The container only accepts traversals in its iteratePairwise() method that implement this interface.
 */
template <class Particle>
class VLCTraversalInterface {
 public:
  /// Verlet list storage
  using verlet_storage_type = std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>;

  /**
   * Sets the verlet list for the traversal to iterate over.
   * @param verlet The verlet list to iterate over.
   */
  virtual void setVerletList(verlet_storage_type &verlet) { _verletList = &verlet; }

 protected:
  /**
   * Iterate over the verlet list of a given cell.
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param verlet
   * @param cellIndex
   * @param pairwiseFunctor
   */
  template <class PairwiseFunctor, bool useNewton3>
  void iterateVerletListsCell(verlet_storage_type &verlet, unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor) {
    for (auto &list : verlet[cellIndex]) {
      Particle &i = *list.first;
      for (auto j_ptr : list.second) {
        Particle &j = *j_ptr;
        pairwiseFunctor->AoSFunctor(i, j, useNewton3);
      }
    }
  }

 protected:
  /**
   * The verlet list to iterate over.
   */
  verlet_storage_type *_verletList;
};

}  // namespace autopas
