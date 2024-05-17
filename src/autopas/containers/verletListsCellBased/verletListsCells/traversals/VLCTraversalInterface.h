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
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 */
template <class Particle, class NeighborList>
class VLCTraversalInterface : public PairwiseTraversalInterface {
 public:
  VLCTraversalInterface() = delete;

  /**
   * Constructor of the VLCTraversalInterface.
   * @param typeOfList indicates the type of neighbor list as an enum value, currently only used for getTraversalType
   */
  VLCTraversalInterface(ContainerOption typeOfList) : _typeOfList(typeOfList) {}

  /**
   * Sets the verlet list for the traversal to iterate over.
   * @param verlet The verlet list to iterate over.
   */
  virtual void setVerletList(NeighborList &verlet) { _verletList = &verlet; }

 protected:
  /**
   * Iterate over all neighbor lists list of a given cell.
   * @tparam PairwiseFunctor
   * @param neighborLists A suitable neighbor list.
   * @param cellIndex
   * @param pairwiseFunctor
   * @param dataLayout
   * @param useNewton3
   */
  template <class PairwiseFunctor>
  void processCellLists(NeighborList &neighborLists, unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor,
                        DataLayoutOption dataLayout, bool useNewton3) {
    processCellListsImpl<PairwiseFunctor>(neighborLists, cellIndex, pairwiseFunctor, dataLayout, useNewton3);
  }

  /**
   * Load the SoA from the respective neighbor list.
   * @tparam pairwiseFunctor
   * @param pairwiseFunctor
   * @param neighborLists
   */
  template <class PairwiseFunctor>
  void loadSoA(PairwiseFunctor *pairwiseFunctor, NeighborList &neighborLists) {
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
  void extractSoA(PairwiseFunctor *pairwiseFunctor, NeighborList &neighborLists) {
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

  /**
   * The type of neighbor list as an enum value.
   */
  ContainerOption _typeOfList;

 private:
  /**
   * Processing of the VLCAllCellsNeighborList type of neighbor list (neighbor list for every cell).
   * @tparam PairwiseFunctor
   * @param neighborList
   * @param cellIndex
   * @param pairwiseFunctor
   * @param dataLayout
   * @param useNewton3
   */
  template <class PairwiseFunctor>
  void processCellListsImpl(VLCAllCellsNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3) {
    if (dataLayout == DataLayoutOption::aos) {
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
        if (not neighbors.empty()) {
          pairwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, neighbors, useNewton3);
        }
      }
    }
  }

  /**
   * Processing of the pairwise Verlet type of neighbor list (neighbor list for every pair of neighboring cells).
   * @tparam PairwiseFunctor
   * @param neighborList
   * @param cellIndex
   * @param pairwiseFunctor
   * @param dataLayout
   * @param useNewton3
   */
  template <class PairwiseFunctor>
  void processCellListsImpl(VLCCellPairNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3) {
    if (dataLayout == DataLayoutOption::aos) {
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

    else if (dataLayout == DataLayoutOption::soa) {
      auto &_soaList = neighborList.getSoANeighborList();
      // iterate over soa and call soaFunctorVerlet for each of the neighbor lists
      for (auto &cellPair : _soaList[cellIndex]) {
        for (auto &[particleIndex, neighbors] : cellPair) {
          if (not neighbors.empty()) {
            pairwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, neighbors, useNewton3);
          }
        }
      }
    }
  }
};

}  // namespace autopas
