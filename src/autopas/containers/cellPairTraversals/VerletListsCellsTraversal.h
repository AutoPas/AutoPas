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
template <class PairwiseFunctor, bool useNewton3>
class VerletListsCellsTraversal {
 public:
  /**
   * Constructor of the verlet traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  VerletListsCellsTraversal(PairwiseFunctor *pairwiseFunctor) : _pairwiseFunctor(pairwiseFunctor) {}

 protected:
  /**
   * Iterate over the verlet list of a given cell.
   * @param verlet
   * @param cellIndex
   */
  template <class Particle>
  inline void iterateVerletListsCell(std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>> &verlet,
                                     unsigned long cellIndex) {
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
