/**
 * @file VerletListsCellsTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"

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
class VerletListsCellsTraversal {
 public:
  /**
   * Shorthand for VerletListsCellsHelpers<Particle>::NeighborListsType.
   */
  using NeighborListsType = typename VerletListsCellsHelpers<Particle>::NeighborListsType;

  /**
   * Sets the verlet list for the traversal to iterate over.
   * @param verlet The verlet list to iterate over.
   */
  virtual void setVerletList(NeighborListsType &verlet) { _verletList = &verlet; }

 protected:
  /**
   * Iterate over all neighbor lists list of a given cell.
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param neighborLists Vector of neighbor lists. One for each particle in the cell.
   * @param cellIndex
   * @param pairwiseFunctor
   */
  template <class PairwiseFunctor, bool useNewton3>
  void processCellLists(NeighborListsType &neighborLists, unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor) {
    for (auto &[particlePtr, neighbors] : neighborLists[cellIndex]) {
      Particle &particle = *particlePtr;
      for (auto neighborPtr : neighbors) {
        Particle &neighbor = *neighborPtr;
        pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
      }
    }
  }

  /**
   * The verlet list to iterate over.
   */
  NeighborListsType *_verletList;
};

}  // namespace autopas
