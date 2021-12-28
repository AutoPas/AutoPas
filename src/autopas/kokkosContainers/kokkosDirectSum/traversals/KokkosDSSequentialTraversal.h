/**
 * @file KokkosDSSequentialTraversal.h
 * @author lgaertner
 * @date 10.11.21
 */

#pragma once

#include <vector>

#include "KokkosDSTraversalInterface.h"
#include "autopas/kokkosContainers/KokkosCellPairTraversals/KokkosCellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/KokkosCellFunctor.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This sum defines the traversal typically used by the DirectSum container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class KokkosDSSequentialTraversal : public KokkosCellPairTraversal<ParticleCell>,
                                    public KokkosDSTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor for the DirectSum traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff (this is enough for the directsum traversal, please don't use the interaction length here.)
   */
  explicit KokkosDSSequentialTraversal(PairwiseFunctor *pairwiseFunctor, double cutoff)
      : KokkosCellPairTraversal<ParticleCell>({2, 1, 1}),
        _cellFunctor(pairwiseFunctor),
        _dataLayoutConverter(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::kokkos_ds_sequential; }

  [[nodiscard]] bool isApplicable() const override { return true; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

  void initTraversal() override {
    // does nothing until aos / soa is supported
  }

  void endTraversal() override {
    // does nothing until aos / soa is supported
  }

  /**
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::KokkosCellFunctor<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, true> _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void KokkosDSSequentialTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
    _cellFunctor.processCell(this->cells(0));
    _cellFunctor.processCellPair(this->cells(0), this->cells(1));
}

}  // namespace autopas
