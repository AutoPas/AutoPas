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
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C08BasedTraversal : public CBasedTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit C08BasedTraversal(const std::array<unsigned long, 3>& dims, PairwiseFunctor* pairwiseFunctor,
                             const double cutoff = 1.0, const std::array<double, 3>& cellLength = {1.0, 1.0, 1.0})
      : CBasedTraversal<ParticleCell>(dims, cutoff, cellLength) {}

  /**
   * C08 traversals are always usable.
   * @return
   */
  bool isApplicable() override { return true; }

 protected:
  /**
   * The main traversal of the C08Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c08Traversal(LoopBody&& loopBody);
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
template <typename LoopBody>
inline void C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::c08Traversal(LoopBody&& loopBody) {
  const auto end = ArrayMath::sub(this->_cellsPerDimension, this->_overlap);
  const auto stride = ArrayMath::addScalar(this->_overlap, 1ul);
  this->cTraversal(std::forward<LoopBody>(loopBody), end, stride);
}
}  // namespace autopas
