/**
 * @file LCSlicedC02Traversal.h
 *
 * @date 24 May 2018
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/SlicedC02BasedTraversal.h"
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
 */
template <class ParticleCell, class PairwiseFunctor>
class LCSlicedC02Traversal : public SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor>,
                             public LCTraversalInterface {
 public:
  /**
   * Constructor of the colored sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit LCSlicedC02Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                double interactionLength, const std::array<double, 3> &cellLength,
                                DataLayoutOption dataLayout, bool useNewton3)
      : SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                               dataLayout, useNewton3, true),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap,
                     dataLayout, useNewton3) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_sliced_c02; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellHandler.setSortingThreshold(sortingThreshold); }

 private:
  LCC08CellHandler<ParticleCell, PairwiseFunctor> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor>
inline void LCSlicedC02Traversal<ParticleCell, PairwiseFunctor>::traverseParticles() {
  auto &cells = *(this->_cells);
  this->cSlicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    auto id = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, id, this->_cellsPerDimension);
  });
}

}  // namespace autopas
