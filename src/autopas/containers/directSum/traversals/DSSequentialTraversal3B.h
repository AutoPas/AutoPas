/**
 * @file DSSequentialTraversal3B.h
 * @author M. Muehlhaeusser
 * @date 07/25/23
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
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class TriwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class DSSequentialTraversal3B : public CellTraversal<ParticleCell>,
                                public DSTraversalInterface<ParticleCell>,
                                public TraversalInterface<InteractionTypeOption::threeBody> {
 public:
  /**
   * Constructor for the DirectSum traversal.
   * @param triwiseFunctor The functor that defines the interaction of three particles.
   * @param cutoff cutoff (this is enough for the directsum traversal, please don't use the interaction length here.)
   */
  explicit DSSequentialTraversal3B(TriwiseFunctor *triwiseFunctor, double cutoff)
      : CellTraversal<ParticleCell>({2, 1, 1}),
        _cellFunctor(triwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/),
        _dataLayoutConverter(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ds_sequential_3b; }

  [[nodiscard]] bool isApplicable() const override { return true; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

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
   * @copydoc autopas::CellTraversal::setUseSorting()
   */
  void setUseSorting(bool useSorting) override { _cellFunctor.setUseSorting(useSorting); }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor3B<typename ParticleCell::ParticleType, ParticleCell, TriwiseFunctor, dataLayout, useNewton3,
                          true>
      _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<TriwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class TriwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void DSSequentialTraversal3B<ParticleCell, TriwiseFunctor, dataLayout, useNewton3>::traverseParticleTriplets() {
  auto &cells = *(this->_cells);
  // Assume cell[0] is the main domain and cell[1] is the halo
  _cellFunctor.processCell(cells[0]);
  _cellFunctor.processCellPair(cells[0], cells[1]);
}

}  // namespace autopas
