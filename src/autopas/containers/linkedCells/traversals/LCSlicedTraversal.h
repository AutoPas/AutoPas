/**
 * @file LCSlicedTraversal.h
 *
 * @date 20 Apr 2018
 * @author gratl
 */

#pragma once

#include <algorithm>

#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/SlicedLockBasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the (locked) sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts the domain into multiple slices along
 * this dimension. Slices are assigned to the threads in a round robin fashion. Each thread locks the cells on the
 * boundary wall to the previous slice with one lock. This lock is lifted as soon the boundary wall is fully processed.
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam Functor_T The functor that defines the interaction of two or three particles.
 */
template <class ParticleCell_T, class Functor_T>
class LCSlicedTraversal : public SlicedLockBasedTraversal<ParticleCell_T, Functor_T>, public LCTraversalInterface {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e., the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction of two or three particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit LCSlicedTraversal(const std::array<unsigned long, 3> &dims, Functor_T *functor, double interactionLength,
                             const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : SlicedLockBasedTraversal<ParticleCell_T, Functor_T>(dims, functor, interactionLength, cellLength, dataLayout,
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
  LCC08CellHandler<ParticleCell_T, Functor_T> _cellHandler;
};

template <class ParticleCell_T, class Functor_T>
inline void LCSlicedTraversal<ParticleCell_T, Functor_T>::traverseParticles() {
  auto &cells = *(this->_cells);
  this->slicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    auto id = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, id);
  });
}

}  // namespace autopas
