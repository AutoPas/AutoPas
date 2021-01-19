/**
 * @file VLCTraversalInterface.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCAllCellsNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"

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
   * @param dataLayout
   */
  template <class PairwiseFunctor, bool useNewton3>
  void processCellLists(NeighborList &neighborLists, unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor,
                        DataLayoutOption::Value dataLayout) {
    processCellListsImpl<PairwiseFunctor, useNewton3>(neighborLists, cellIndex, pairwiseFunctor, dataLayout);
  }

  /**
   * Load the SoA from the respective neighbor list.
   * @tparam pairwiseFunctor
   * @param pairwiseFunctor
   * @param neighborLists
   */
  template <class PairwiseFunctor>
  void setupLoadSoA(PairwiseFunctor *pairwiseFunctor, NeighborList &neighborLists) {
    // send to loadSoA in the neighbor list
    _soa = neighborLists.loadSoA(pairwiseFunctor);
  }

  /**
   * Extract the SoA from the respective neighbor list.
   * @tparam pairwiseFunctor
   * @param pairwiseFunctor
   * @param neighborLists
   */
  template <class PairwiseFunctor>
  void setupExtractSoA(PairwiseFunctor *pairwiseFunctor, NeighborList &neighborLists) {
    // send to extractSoA in the neighbor list
    neighborLists.extractSoA(pairwiseFunctor);
    _soa = nullptr;
  }

  /**
   * The verlet list to iterate over.
   */
  NeighborList *_verletList;

  /**
   * Structure of arrays to be used if the data layout is SoA.
   */
  SoA<typename Particle::SoAArraysType> *_soa;

 private:
  /** Processing of the VLCAllCellsNeighborList type of neighbor list (neighbor list for every cell).
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param neighborList
   * @param cellIndex
   * @param pairwiseFunctor
   * @param dataLayout
   */
  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(VLCAllCellsNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value dataLayout) {
    if (dataLayout == DataLayoutOption::Value::aos) {
      auto &internalList = neighborList.getAoSNeighborList();
      for (auto &[particlePtr, neighbors] : internalList[cellIndex]) {
        Particle &particle = *particlePtr;
        for (auto neighborPtr : neighbors) {
          Particle &neighbor = *neighborPtr;
          pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
        }
      }
    }

    else if (dataLayout == DataLayoutOption::soa) {
      auto &_soaList = neighborList.getSoANeighborList();
      for (auto &[particleIndex, neighbors] : _soaList[cellIndex]) {
        pairwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, neighbors, useNewton3);
      }
    }
  }

  /** Processing of the pairwise Verlet type of neighbor list (neighbor list for every pair of neighboring cells).
   * @tparam PairwiseFunctor
   * @tparam useNewton3
   * @param neighborList
   * @param cellIndex
   * @param pairwiseFunctor
   * @param dataLayout
   */
  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(VLCCellPairNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value dataLayout) {
    if (dataLayout == DataLayoutOption::Value::aos) {
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

    else if (dataLayout == DataLayoutOption::Value::soa) {
      auto &_soaList = neighborList.getSoANeighborList();
      // iterate over soa and call soaFunctorVerlet for each of the neighbor lists
      for (auto &cellPair : _soaList[cellIndex]) {
        for (auto &[particleIndex, neighbors] : cellPair) {
          if(!neighbors.empty())
          {
            pairwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, neighbors, useNewton3);
          }
        }
      }
    }
  }
};

}  // namespace autopas
