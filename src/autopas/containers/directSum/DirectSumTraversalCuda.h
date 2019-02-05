/**
 * @file DirectSumTraversalCuda.h
 * @author jspahl
 * @date 12/9/18
 */

#pragma once

#include <vector>
#include "DirectSumTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctorCuda.h"

namespace autopas {
/**
 * This sum defines the traversal typically used by the DirectSum container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class DirectSumTraversalCuda : public CellPairTraversal<ParticleCell>, public DirectSumTraversalInterface<ParticleCell> {

 public:
  /**
   * Constructor for the DirectSum traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  DirectSumTraversalCuda(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        _cellFunctor(
            CellFunctorCuda<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3>(
                pairwiseFunctor)) {}

  TraversalOptions getTraversalType() override;

  bool isApplicable() override;

  /**
   * @copydoc CellPairTraversal::traverseCellPairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctorCuda<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3> _cellFunctor;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
TraversalOptions DirectSumTraversalCuda<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::getTraversalType() {
  return TraversalOptions::directSumTraversalCuda;
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
bool DirectSumTraversalCuda<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::isApplicable() {
  return true;
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
void DirectSumTraversalCuda<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  // Assume cell[0] is the main domain and cell[1] is the halo
  _cellFunctor.processCell(cells[0]);
  _cellFunctor.processCellPair(cells[0], cells[1]);
}
}  // namespace autopas
