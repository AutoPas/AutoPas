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

template <class Particle>
class VerletListsCellsNeighborList;

template <class Particle>
class PairwiseVerletNeighborList;

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
  virtual void setVerletList(NeighborList &verlet) {
    _verletList = &verlet;
    internalList = _verletList->getAoSNeighborList();  // only this->verletList is being passed to processCellLists so
                                                       // internalList should be safe to use below
    _soaList = _verletList->getSoANeighborList();
  }

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
  void processCellLists(NeighborList &neighborLists, unsigned long cellIndex, PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value dataLayout) {
    processCellListsImpl<PairwiseFunctor, useNewton3>(neighborLists, cellIndex, pairwiseFunctor, dataLayout);
  }

  //TODO
  template <class PairwiseFunctor>
  void setupLoadSoA(PairwiseFunctor *pairwiseFunctor)
  {
    //send to loadSoA in PairwiseVerletNeighborList
    _soa = _verletList->loadSoA(pairwiseFunctor);
  }

  //TODO
  template <class PairwiseFunctor>
  void setupExtractSoA(PairwiseFunctor *pairwiseFunctor)
  {
    //send to extractSoA in PaiwiseVL
    _verletList->extractSoA(pairwiseFunctor);
    _soa = nullptr;
  }

  /**
   * The neighbor list class which contains the verlet list to iterate over.
   */
  NeighborList *_verletList;
  /**
   * The verlet list to iterate over.
   */
  typename NeighborList::listType internalList;
  typename NeighborList::soaListType _soaList;

  //TODO
  SoA<typename Particle::SoAArraysType>* _soa;

 private:
  /** Processing of the VerletListsCellsNeighborList type of neighbor list (neighbor list for every cell).*/
  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(VerletListsCellsNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value dataLayout) {

    //if(dataLayout == DataLayoutOption::Value::aos)
    //{
      for (auto &[particlePtr, neighbors] : internalList[cellIndex]) {
        Particle &particle = *particlePtr;
        for (auto neighborPtr : neighbors) {
          Particle &neighbor = *neighborPtr;
          pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
        }
      }
    //}
  }

  /** Processing of the pairwise Verlet type of neighbor list (neighbor list for every pair of neighboring cells).*/
  template <class PairwiseFunctor, bool useNewton3>
  void processCellListsImpl(PairwiseVerletNeighborList<Particle> &neighborList, unsigned long cellIndex,
                            PairwiseFunctor *pairwiseFunctor, DataLayoutOption::Value dataLayout) {
    /*if(dataLayout == DataLayoutOption::Value::aos)
    {
      size_t numberOfCellToInteract = 27;
      for (size_t neighborCellIndex = 0; neighborCellIndex < numberOfCellToInteract; neighborCellIndex++) {
        for (auto &[particlePtr, neighbors] : internalList[cellIndex][neighborCellIndex]) {
          Particle &particle = *particlePtr;
          for (auto neighborPtr : neighbors) {
            Particle &neighbor = *neighborPtr;
            pairwiseFunctor->AoSFunctor(particle, neighbor, useNewton3);
          }
        }
      }
    }*/

    //else if (dataLayout == DataLayoutOption::Value::soa)
    //{
      //iterate over soa and call soaFunctorVerlet for  each of the small neighbor lists
      for(size_t secondCell = 0; secondCell < _soaList[cellIndex].size(); secondCell++)
      {
        for(size_t particleIndex = 0; particleIndex < _soaList[cellIndex][secondCell].size(); particleIndex++)
        {
          std::cout << "called " <<std::endl;
          pairwiseFunctor->SoAFunctorVerlet(*_soa, _soaList[cellIndex][secondCell][particleIndex].first, _soaList[cellIndex][secondCell][particleIndex].second, useNewton3);
        }
      }
    //}
  }
};

}  // namespace autopas
