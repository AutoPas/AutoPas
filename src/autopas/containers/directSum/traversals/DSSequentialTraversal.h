/**
 * @file DSSequentialTraversal.h
 * @author F. Gratl
 * @date 11/23/18
 */

#pragma once

#include <vector>

#include "DSTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This sum defines the traversal typically used by the DirectSum container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class DSSequentialTraversal : public CellPairTraversal<ParticleCell>, public DSTraversalInterface<ParticleCell> {
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
      : CellPairTraversal<ParticleCell>({2, 1, 1}, dataLayout, useNewton3),
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
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override;

  /**
   * @copydoc autopas::CellPairTraversal::setSortingThreshold()
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
  //  Process interactions with all six halo cells
  _cellFunctor.processCellPair(cells[0], cells[1], std::array<double, 3>({-1., 0., 0.}));
  _cellFunctor.processCellPair(cells[0], cells[2], std::array<double, 3>({1., 0., 0.}));
  _cellFunctor.processCellPair(cells[0], cells[3], std::array<double, 3>({0., -1., 0.}));
  _cellFunctor.processCellPair(cells[0], cells[4], std::array<double, 3>({0., 1., 0.}));
  _cellFunctor.processCellPair(cells[0], cells[5], std::array<double, 3>({0., 0., -1.}));
  _cellFunctor.processCellPair(cells[0], cells[6], std::array<double, 3>({0., 0., 1.}));
}

}  // namespace autopas
