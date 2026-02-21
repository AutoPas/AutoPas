/**
 * @file VLCTraversalInterface.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

#include "autopas/utils/checkFunctorType.h"
#include "autopas/options/ContainerOption.h"

namespace autopas {

/**
 *
 * Forward declaration of the neighbor list types used in this interface.
 * This avoids circular dependencies in the header files.
 * VLCAllCellsNeighborList.h/VLCCellPairNeighborList.h -> VLCTraversalInterface.h (This file)
 */
template <class Particle_T>
class VLCAllCellsNeighborList;

template <class Particle_T>
class VLCCellPairNeighborList;

/**
 * This class provides the Traversal Interface for the verlet lists cells container.
 *
 * This class handles traversals through the cell structures with neighbor lists.
 * Derived classes handle the order through which the cells are traversed.
 *
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 */
template <class Particle_T, class NeighborList>
class VLCTraversalInterface {
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
  template <class Functor>
  void processCellLists(NeighborList &neighborLists, unsigned long cellIndex, Functor *functor,
                        DataLayoutOption dataLayout, bool useNewton3) {
    if constexpr (utils::isPairwiseFunctor<Functor>()) {
      processCellListsPairwiseImpl<Functor>(neighborLists, cellIndex, functor, dataLayout, useNewton3);
    } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
      processCellListsTriwiseImpl<Functor>(neighborLists, cellIndex, functor, dataLayout, useNewton3);
    } else {
      utils::ExceptionHandler::exception(
          "VLCTraversalInterface::processCellLists(): Functor {} is not of type PairwiseFunctor or TriwiseFunctor.",
          functor->getName());
    }
  }

  /**
   * Load the SoA from the respective neighbor list.
   * @tparam PairwiseFunctor
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
   * @tparam PairwiseFunctor
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
  SoA<typename Particle_T::SoAArraysType> *_soa;

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
  void processCellListsPairwiseImpl(VLCAllCellsNeighborList<Particle_T> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3) {
    if (dataLayout == DataLayoutOption::aos) {
      auto &aosList = neighborList.getAoSNeighborList();
      for (auto &[particlePtr, neighbors] : aosList[cellIndex]) {
        Particle_T &particle = *particlePtr;
        for (auto neighborPtr : neighbors) {
          Particle_T &neighbor = *neighborPtr;
          pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
        }
      }
    }

    else if (dataLayout == DataLayoutOption::soa) {
      auto &soaList = neighborList.getSoANeighborList();
      for (auto &[particleIndex, neighbors] : soaList[cellIndex]) {
        if (not neighbors.empty()) {
          pairwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, neighbors, useNewton3);
        }
      }
    }
  }

  template <class TriwiseFunctor>
void processCellListsTriwiseImpl(VLCAllCellsNeighborList<Particle_T> &neighborList, unsigned long cellIndex,
                          TriwiseFunctor *triwiseFunctor, DataLayoutOption dataLayout, bool useNewton3) {
    if (dataLayout == DataLayoutOption::aos) {
      auto &aosList = neighborList.getAoSNeighborList();
      for (auto &[particlePtr, neighbors] : aosList[cellIndex]) {
        Particle_T &particle = *particlePtr;
        for (auto neighborPtr1 = neighbors.begin(); neighborPtr1 != neighbors.end(); ++neighborPtr1) {
          Particle_T &neighbor1 = **neighborPtr1;
            for (auto neighborPtr2 = std::next(neighborPtr1); neighborPtr2 != neighbors.end(); ++neighborPtr2) {
              Particle_T &neighbor2 = **neighborPtr2;
              triwiseFunctor->AoSFunctor(particle, neighbor1, neighbor2, useNewton3);
            }
        }
      }
    }

    else if (dataLayout == DataLayoutOption::soa) {
      auto &soaList = neighborList.getSoANeighborList();
      for (auto &[particleIndex, neighbors] : soaList[cellIndex]) {
        if (not neighbors.empty()) {
          triwiseFunctor->SoAFunctorVerlet(*_soa, particleIndex, neighbors, useNewton3);
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
  void processCellListsPairwiseImpl(VLCCellPairNeighborList<Particle_T> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3) {
    if (dataLayout == DataLayoutOption::aos) {
      auto &aosList = neighborList.getAoSNeighborList();
      for (auto &cellPair : aosList[cellIndex]) {
        for (auto &[particlePtr, neighbors] : cellPair) {
          Particle_T &particle = *particlePtr;
          for (auto neighborPtr : neighbors) {
            Particle_T &neighbor = *neighborPtr;
            pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
          }
        }
      }
    }

    else if (dataLayout == DataLayoutOption::soa) {
      auto &soaList = neighborList.getSoANeighborList();
      // iterate over soa and call soaFunctorVerlet for each of the neighbor lists
      for (auto &cellPair : soaList[cellIndex]) {
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
