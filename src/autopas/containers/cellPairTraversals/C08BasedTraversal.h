/**
 * @file C08BasedTraversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CBasedTraversal.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The traversal is defined in the function c08Traversal and uses 8 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID in each direction are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class C08BasedTraversal : public CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cell block, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit C08BasedTraversal(const std::array<uint64_t, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                             const double interactionLength, const std::array<double, 3> &cellLength)
      : CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor, interactionLength,
                                                                               cellLength) {}

 protected:
  /**
   * The main traversal of the C08Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c08Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <typename LoopBody>
inline void C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::c08Traversal(
    LoopBody &&loopBody) {
  const auto end = utils::ArrayMath::sub(this->_cellsPerDimension, this->_overlap);
  const auto stride = utils::ArrayMath::addScalar(this->_overlap, 1ul);
  this->cTraversal(std::forward<LoopBody>(loopBody), end, stride);
}
}  // namespace autopas
