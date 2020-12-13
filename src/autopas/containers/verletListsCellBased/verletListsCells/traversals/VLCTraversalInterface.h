/**
 * @file VLCTraversalInterface.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

#include "autopas/containers/verletListsCellBased/verletListsCells/PairwiseVerletNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsNeighborList.h"

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
   * Sets the verlet list for the traversal to iterate over.
   * @param verlet The verlet list to iterate over.
   */
  virtual void setVerletList(NeighborList &verlet) { _verletList = &verlet; }

 protected:
  /**
   * Iterate over all neighbor lists list of a given cell.
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param neighborLists A suitable neighbor list.
   * @param cellIndex
   * @param pairwiseFunctor
   */
  template <class PairwiseFunctor, bool useNewton3>
  void processCellLists(NeighborList &neighborLists, unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor) {
    processCellListsImpl<PairwiseFunctor, useNewton3>(neighborLists, cellIndex, pairwiseFunctor);
  }

  /**
   * The verlet list to iterate over.
   */
  NeighborList *_verletList;

 private:
  /** Processing of the VerletListsCellsNeighborList type of neighbor list (neighbor list for every cell).
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param neighborList
   * @param cellIndex
   * @param pairwiseFunctor
   */
  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(VerletListsCellsNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor) {
    auto &internalList = neighborList.getAoSNeighborList();
    for (auto &[particlePtr, neighbors] : internalList[cellIndex]) {
      Particle &particle = *particlePtr;
      for (auto neighborPtr : neighbors) {
        Particle &neighbor = *neighborPtr;
        pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
      }
    }
  }

  /** Processing of the pairwise Verlet type of neighbor list (neighbor list for every pair of neighboring cells).
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param neighborList
   * @param cellIndex
   * @param pairwiseFunctor
   */
  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(PairwiseVerletNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor) {
    auto &internalList = neighborList.getAoSNeighborList();
    for (auto &cellPair : internalList[cellIndex]) {
      for (auto &[particlePtr, neighbors] : cellPair) {
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
