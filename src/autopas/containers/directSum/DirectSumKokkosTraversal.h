/**
 * @file DirectSumKokkosTraversal.h
 * @author M. Geitner
 * @date 04/18/2019
 */

#pragma once

#include <vector>
#include "DirectSumTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"

#ifdef AUTOPAS_KOKKOS
#include <Kokkos_Core.hpp>
#endif

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
class DirectSumKokkosTraversal : public CellPairTraversal<ParticleCell>,
                                 public DirectSumTraversalInterface<ParticleCell> {
  /**
   * Constructor for the DirectSum traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
 public:
  DirectSumKokkosTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        _cellFunctor(
            CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3>(
                pairwiseFunctor)) {}

  TraversalOption getTraversalType() override { return TraversalOption::directSumKokkosTraversal; }

  bool isApplicable() override { return true; }
  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3> _cellFunctor;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
void DirectSumKokkosTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  // Assume cell[0] is the main domain and cell[1] is the halo

#ifdef AUTOPAS_KOKKOS
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i){

                              //@FIXME Kokkos cannot use non-const functions
                              // typedef Kokkos::View<CellFunctor<typename ParticleCell::ParticleType, ParticleCell,
                              // PairwiseFunctor, useSoA, useNewton3>(
                              //              pairwiseFunctor)*>   ViewVectorType;
                              // ViewVectorType functor( "functor", 1);
                              //_cellFunctor.processCell(cells[0]);
                              //_cellFunctor.processCellPair(cells[0], cells[1]);

                          });
#endif
}

}  // namespace autopas