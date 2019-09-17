/**
 * @file C18BasedTraversal.h
 * @author nguyen
 * @date 06.09.18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CBasedTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c18 base step.
 * The traversal is defined in the function c18Traversal and uses 18 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class C18BasedTraversal : public CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit C18BasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                             const double interactionLength, const std::array<double, 3> &cellLength)
      : CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor, interactionLength,
                                                                               cellLength) {}

 protected:
  /**
   * The main traversal of the C18Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   * @tparam allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * it will not iterate over all cells and instead skip the last few cells, because they will be covered by the base
   * step.
   */
  template <bool allCells = false, typename LoopBody>
  inline void c18Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
template <bool allCells, typename LoopBody>
inline void C18BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::c18Traversal(
    LoopBody &&loopBody) {
  const std::array<unsigned long, 3> stride = {2ul * this->_overlap[0] + 1ul, 2ul * this->_overlap[1] + 1ul,
                                               this->_overlap[2] + 1ul};
  auto end(this->_cellsPerDimension);
  if (not allCells) {
    end[2] -= this->_overlap[2];
  }
  this->cTraversal(std::forward<LoopBody>(loopBody), end, stride);
}

}  // namespace autopas
