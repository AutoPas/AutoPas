/**
 * @file VerletListsCellsTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>
#include "autopas/containers/VerletListsCellsHelpers.h"

namespace autopas {

/**
 * A verlet list traversal.
 * This class handles traversals through the cell structures with neighbor lists.
 * Derived classes handle the order through which the cells are traversed.
 */
template <class Particle, class PairwiseFunctor, bool useNewton3>
class VerletListsCellsTraversal {
 public:
  typedef typename VerletListsCellsHelpers<Particle>::VerletList_storage_type verlet_storage_type;

  /**
   * Constructor of the verlet traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  VerletListsCellsTraversal(PairwiseFunctor *pairwiseFunctor) : _pairwiseFunctor(pairwiseFunctor) {}

  /**
   * Traverse verlet lists of all cells.
   * This function needs to be implemented by derived classes.
   * @param verlet verlet lists for each cell
   */
  virtual void traverseCellVerlet(verlet_storage_type &verlet) = 0;

 protected:
  /**
   * Iterate over the verlet list of a given cell.
   * @param verlet
   * @param cellIndex
   */
  inline void iterateVerletListsCell(verlet_storage_type &verlet, unsigned long cellIndex) {
    for (auto &list : verlet[cellIndex]) {
      Particle &i = *list.first;
      for (auto j_ptr : list.second) {
        Particle &j = *j_ptr;
        _pairwiseFunctor->AoSFunctor(i, j, useNewton3);
      }
    }
  }

 private:
  /**
   * PairwiseFunctor to be used for the traversal defining the interaction between two particles.
   */
  PairwiseFunctor *_pairwiseFunctor;
};

}  // namespace autopas
