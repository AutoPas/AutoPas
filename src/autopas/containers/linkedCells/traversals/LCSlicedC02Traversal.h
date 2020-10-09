/**
 * @file LCSlicedC02Traversal.h
 *
 * @date 24 May 2018
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "LCTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/SlicedC02BasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the colored sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into as many slices as possible along this dimension. Unlike the regular
 * sliced traversal, this version uses a 2-coloring to prevent race conditions, instead of
 * locking the starting layers.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCSlicedC02Traversal : public SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                             public LCTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the colored sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit LCSlicedC02Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                const double interactionLength, const std::array<double, 3> &cellLength)
      : SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(
            dims, pairwiseFunctor, interactionLength, cellLength, true),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

  void traverseParticlePairs() override;

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_sliced_c02; }

 private:
  LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCSlicedC02Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  this->cSlicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    auto id = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, id);
  });
}

}  // namespace autopas
