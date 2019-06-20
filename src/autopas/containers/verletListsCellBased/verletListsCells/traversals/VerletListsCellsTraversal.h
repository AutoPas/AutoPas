/**
 * @file VerletListsCellsTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

namespace autopas {

/**
 * A verlet list traversal.
 * This class handles traversals through the cell structures with neighbor lists.
 * Derived classes handle the order through which the cells are traversed.
 */
template <class Particle>
class VerletListsCellsTraversal {
 public:
  /// Verlet list storage
  typedef std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>> verlet_storage_type;

  /**
   * Sets the verlet list for the traversal to iterate over.
   * @param verlet The verlet list to iterate over.
   */
  virtual void setVerletList(verlet_storage_type &verlet) { _verletList = &verlet; }

 protected:
  /**
   * Iterate over the verlet list of a given cell.
   * @param verlet
   * @param cellIndex
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
  verlet_storage_type *_verletList;
};

}  // namespace autopas
