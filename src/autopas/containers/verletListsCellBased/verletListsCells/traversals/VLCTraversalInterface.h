/**
 * @file VLCTraversalInterface.h
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
template <class Particle, class NeighborList>
class VLCTraversalInterface {
 public:
  /**
   * Shorthand for VerletListsCellsHelpers<Particle>::NeighborListsType.
   */
  using NeighborListsType = typename VerletListsCellsHelpers<Particle>::NeighborListsType;
  /**
   * Shorthand for VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType.
   */
  using PairwiseNeighborListsType = typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType;

  /**
   * Sets the verlet list for the traversal to iterate over.
   * @param verlet The verlet list to iterate over.
   */
  virtual void setVerletList(NeighborList &verlet) { _verletList = &verlet; }

 protected:
  /**
   * Iterate over all neighbor lists list of a given cell.
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param neighborLists Vector of neighbor lists. One for each particle in the cell.
   * @param cellIndex
   * @param pairwiseFunctor
   */
  template <class PairwiseFunctor, bool useNewton3, class NL>
  void processCellLists(NL &neighborLists, unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor) {
    processCellListsImpl<PairwiseFunctor, useNewton3>(neighborLists, cellIndex, pairwiseFunctor);
  }

  /**
   * The verlet list to iterate over.
   */
  NeighborList *_verletList;

 private:
  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(NeighborListsType &neighborLists, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor) {
    for (auto &[particlePtr, neighbors] : neighborLists[cellIndex]) {
      Particle &particle = *particlePtr;
      for (auto neighborPtr : neighbors) {
        Particle &neighbor = *neighborPtr;
        pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
      }
    }
  }

  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(PairwiseNeighborListsType &neighborLists, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor) {
    size_t numberOfCellToInteract = 27;
    for (size_t neighborCellIndex = 0; neighborCellIndex < numberOfCellToInteract; neighborCellIndex++) {
      for (auto &[particlePtr, neighbors] : neighborLists[cellIndex][neighborCellIndex]) {
        Particle &particle = *particlePtr;
        for (auto neighborPtr : neighbors) {
          Particle &neighbor = *neighborPtr;
          pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
        }
      }
    }
  }
};

}  // namespace autopas
