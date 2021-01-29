/**
 * @file LCSlicedBalancedTraversal.h
 *
 * @date 27 Apr 2020
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "LCTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/SlicedBalancedBasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the balanced sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain along this dimension in one slice (block) per thread. The cut
 * positions are calculated to even out load among the threads. Slices are
 * assigned to the threads in a round robin fashion. Each thread locks the cells
 * on the boundary wall to the previous slice with one lock. This lock is lifted
 * as soon the boundary wall is fully processed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCSlicedBalancedTraversal
    : public SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, true>,
      public LCTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the balanced sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit LCSlicedBalancedTraversal(const std::array<uint64_t, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                     const double interactionLength, const std::array<double, 3> &cellLength)
      : SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, true>(
            dims, pairwiseFunctor, interactionLength, cellLength),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

  void traverseParticlePairs() override;

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_sliced_balanced; }

 private:
  LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCSlicedBalancedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  this->slicedTraversal([&](uint64_t x, uint64_t y, uint64_t z) {
    auto id = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, id);
  });
}

}  // namespace autopas
