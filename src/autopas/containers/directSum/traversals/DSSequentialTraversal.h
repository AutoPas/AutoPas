/**
 * @file DSSequentialTraversal.h
 * @author F. Gratl
 * @date 11/23/18
 */

#pragma once

#include <vector>

#include "DSPairTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This sum defines the traversal typically used by the DirectSum container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class DSSequentialTraversal : public CellTraversal<ParticleCell>, public DSPairTraversalInterface {
 public:
  /**
   * Constructor for the DirectSum traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff (this is enough for the directsum traversal, please don't use the interaction length here.)
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit DSSequentialTraversal(PairwiseFunctor *pairwiseFunctor, double cutoff, DataLayoutOption dataLayout,
                                 bool useNewton3)
      : CellTraversal<ParticleCell>({2, 1, 1}),
        TraversalInterface(dataLayout, useNewton3),
        _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/, dataLayout,
                     useNewton3),
        _dataLayoutConverter(pairwiseFunctor, dataLayout) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ds_sequential; }

  [[nodiscard]] bool isApplicable() const override { return true; }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    for (auto &cell : cells) {
      _dataLayoutConverter.loadDataLayout(cell);
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    for (auto &cell : cells) {
      _dataLayoutConverter.storeDataLayout(cell);
    }
  }

  /**
   * @copydoc PairwiseTraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override;

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellFunctor.setSortingThreshold(sortingThreshold); }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<ParticleCell, PairwiseFunctor, /*bidirectional*/ true> _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor>
void DSSequentialTraversal<ParticleCell, PairwiseFunctor>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  // Assume cell[0] is the main domain and cell[1] is the halo
  _cellFunctor.processCell(cells[0]);
  _cellFunctor.processCellPair(cells[0], cells[1]);
}

}  // namespace autopas
