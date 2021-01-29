/**
 * @file C04BasedTraversal.h
 * @author C. Menges
 * @date 02.06.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CBasedTraversal.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c04 base step.
 *
 * The traversal is defined in the function c04Traversal and uses 4 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID in each direction are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          int collapseDepth = 3>
class C04BasedTraversal : public CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, collapseDepth> {
 public:
  /**
   * Constructor of the c04 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length.
   * @param cellLength cell length.
   */
  explicit C04BasedTraversal(const std::array<uint64_t, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                             const double interactionLength, const std::array<double, 3> &cellLength)
      : CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, collapseDepth>(
            dims, pairwiseFunctor, interactionLength, cellLength) {}

 protected:
  /**
   * The main traversal of the C04Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c04Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          int collapseDepth>
template <typename LoopBody>
inline void C04BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, collapseDepth>::c04Traversal(
    LoopBody &&loopBody) {
  const auto end = utils::ArrayMath::sub(this->_cellsPerDimension, this->_overlap);
  auto stride = utils::ArrayMath::addScalar(this->_overlap, 1ul);
  stride[0] = 1;
  this->cTraversal(std::forward<LoopBody>(loopBody), end, stride);
}
}  // namespace autopas
