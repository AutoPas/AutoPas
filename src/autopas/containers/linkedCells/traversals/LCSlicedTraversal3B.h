/**
 * @file LCSlicedTraversal3B.h
 *
 * @date 21 Nov 2023
 * @author N. Deng
 */

#pragma once

#include <algorithm>

#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/SlicedLockBasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the (locked) sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into multiple slices along this dimension. Slices are
 * assigned to the threads in a round robin fashion. Each thread locks the cells
 * on the boundary wall to the previous slice with one lock. This lock is lifted
 * as soon the boundary wall is fully processed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of three particles.
 */
template <class ParticleCell, class Functor>
class LCSlicedTraversal3B : public SlicedLockBasedTraversal<ParticleCell, Functor>, public LCTraversalInterface {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction of three particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit LCSlicedTraversal3B(const std::array<unsigned long, 3> &dims, Functor *functor,
                               const double interactionLength, const std::array<double, 3> &cellLength,
                               DataLayoutOption dataLayout, bool useNewton3)
      : SlicedLockBasedTraversal<ParticleCell, Functor>(dims, functor, interactionLength, cellLength, dataLayout,
                                                        useNewton3, true),
        _cellHandler(functor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap, dataLayout,
                     useNewton3) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_sliced; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellHandler.setSortingThreshold(sortingThreshold); }

 private:
  LCC08CellHandler<ParticleCell, Functor> _cellHandler;
};

template <class ParticleCell, class Functor>
inline void LCSlicedTraversal3B<ParticleCell, Functor>::traverseParticles() {
  auto &cells = *(this->_cells);
  this->slicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, baseIndex);
  });
}

}  // namespace autopas
