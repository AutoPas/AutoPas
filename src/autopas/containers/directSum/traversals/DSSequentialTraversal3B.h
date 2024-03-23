/**
 * @file DSSequentialTraversal3B.h
 * @author M. Muehlhaeusser
 * @date 25.07.2023
 */

#pragma once

#include <vector>

#include "DSTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This sum defines the traversal typically used by the DirectSum container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam TriwiseFunctor The functor that defines the interaction of three particles.

 */
template <class ParticleCell, class TriwiseFunctor>
class DSSequentialTraversal3B : public CellTraversal<ParticleCell>,
                                public DSTraversalInterface<ParticleCell>,
                                public TraversalInterface<InteractionTypeOption::threeBody> {
 public:
  /**
   * Constructor for the DirectSum traversal.
   * @param triwiseFunctor The functor that defines the interaction of three particles.
   * @param cutoff cutoff (this is enough for the directsum traversal, please don't use the interaction length here.)
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit DSSequentialTraversal3B(TriwiseFunctor *triwiseFunctor, double cutoff, DataLayoutOption dataLayout,
                                   bool useNewton3)
      : CellTraversal<ParticleCell>({2, 1, 1}),
        TraversalInterface<InteractionTypeOption::threeBody>(dataLayout, useNewton3),
        _cellFunctor(triwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/, dataLayout,
                     useNewton3),
        _dataLayoutConverter(triwiseFunctor, dataLayout) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ds_sequential_3b; }

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
   * @copydoc TraversalInterface::traverseParticleTriplets()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticleTriplets() override;

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellFunctor.setSortingThreshold(sortingThreshold); }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor3B<ParticleCell, TriwiseFunctor, true> _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<TriwiseFunctor> _dataLayoutConverter;
};

template <class ParticleCell, class TriwiseFunctor>
void DSSequentialTraversal3B<ParticleCell, TriwiseFunctor>::traverseParticleTriplets() {
  auto &cells = *(this->_cells);
  // Assume cell[0] is the main domain and cell[1] is the halo
  _cellFunctor.processCell(cells[0]);
  _cellFunctor.processCellPair(cells[0], cells[1]);
}

}  // namespace autopas
