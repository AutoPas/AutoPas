/**
 * @file C01BasedTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <autopas/utils/WrapOpenMP.h>
#include "autopas/containers/cellPairTraversals/CBasedTraversal.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c01 base step.
 *
 * The traversal is defined in the function c01Traversal and uses 1 color. Interactions between two cells are allowed
 * only if particles of the first cell are modified. This means that newton3 optimizations are NOT allowed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C01BasedTraversal : public CBasedTraversal<ParticleCell, PairwiseFunctor, DataLayout> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit C01BasedTraversal(const std::array<unsigned long, 3>& dims, PairwiseFunctor* pairwiseFunctor,
                             double cutoff = 1.0, const std::array<double, 3>& cellLength = {1.0, 1.0, 1.0})
      : CBasedTraversal<ParticleCell, PairwiseFunctor, DataLayout>(dims, pairwiseFunctor, cutoff, cellLength) {}

 protected:
  /**
   * The main traversal of the C01Traversal.
   * This provides the structure of the loops and its parallelization.
   * @tparam LoopBody
   * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
   * (x,y,z). If you need additional input from outside, please use captures (by reference).
   */
  template <typename LoopBody>
  inline void c01Traversal(LoopBody&& loopBody);
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
template <typename LoopBody>
inline void C01BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::c01Traversal(
    LoopBody&& loopBody) {
  const auto offset = this->_overlap;
  const auto end = ArrayMath::sub(this->_cellsPerDimension, this->_overlap);
  this->cTraversal(std::forward<LoopBody>(loopBody), end, {1ul, 1ul, 1ul}, offset);
}
}  // namespace autopas
