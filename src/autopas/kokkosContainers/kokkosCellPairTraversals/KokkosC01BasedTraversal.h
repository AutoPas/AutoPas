/**
* @file KokkosC01BasedTraversal.h
* @author lgaertner
* @date 29.12.21
*/

#pragma once

#include "autopas/kokkosContainers/kokkosCellPairTraversals/KokkosCBasedTraversal.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
* This class provides the base for traversals using the c01 base step.
*
* The traversal is defined in the function c01Traversal and uses 1 color. Interactions between two cells are allowed
* only if particles of the first cell are modified. This means that newton3 optimizations are NOT allowed.
*
* @tparam ParticleCell the type of cells
* @tparam PairwiseFunctor The functor that defines the interaction of two particles.
* @tparam dataLayout indicates usage of SoA
*/
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class KokkosC01BasedTraversal : public KokkosCBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit KokkosC01BasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                             double interactionLength, const std::array<double, 3> &cellLength)
      : KokkosCBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(
      dims, pairwiseFunctor, interactionLength, cellLength) {}

 protected:
  /**
   * The main traversal of the C01Traversal.
   * This provides the structure of the loops and its parallelization.
   * @tparam LoopBody
   * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
   * (x,y,z). If you need additional input from outside, please use captures (by reference).
   */
  template <typename LoopBody>
  inline void c01Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <typename LoopBody>
inline void KokkosC01BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::c01Traversal(
    LoopBody &&loopBody) {
  const auto offset = this->_overlap;
  const auto end = utils::ArrayMath::sub(this->_cellsPerDimension, this->_overlap);
  this->cTraversal(std::forward<LoopBody>(loopBody), end, {1ul, 1ul, 1ul}, offset);
}
}  // namespace autopas
