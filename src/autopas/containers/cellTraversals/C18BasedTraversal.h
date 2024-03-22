/**
 * @file C18BasedTraversal.h
 * @author nguyen
 * @date 06.09.18
 */

#pragma once

#include "ColorBasedTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c18 base step.
 * The traversal is defined in the function c18Traversal and uses 18 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType>
class C18BasedTraversal : public ColorBasedTraversal<ParticleCell, Functor, interactionType> {
 public:
  /**
   * Constructor of the lc_c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit C18BasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor, const double interactionLength,
                             const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : ColorBasedTraversal<ParticleCell, Functor, interactionType>(
            dims, functor, interactionLength, cellLength, dataLayout, useNewton3) {}

 protected:
  /**
   * The main traversal of the C18Traversal.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   *
   * @tparam allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the c18 step if allCells is false, iteration will not occur over the last layer of
   * cells (for overlap=1) (in x, y and z direction).
   */
  template <bool allCells = false, typename LoopBody>
  inline void c18Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType>
template <bool allCells, typename LoopBody>
inline void C18BasedTraversal<ParticleCell, Functor, interactionType>::c18Traversal(
    LoopBody &&loopBody) {
  const std::array<unsigned long, 3> stride = {2ul * this->_overlap[0] + 1ul, 2ul * this->_overlap[1] + 1ul,
                                               this->_overlap[2] + 1ul};
  auto end(this->_cellsPerDimension);
  if (not allCells) {
    end[2] -= this->_overlap[2];
  }
  this->colorTraversal(std::forward<LoopBody>(loopBody), end, stride);
}

}  // namespace autopas
