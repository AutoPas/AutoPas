/**
 * @file C08BasedTraversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "ColorBasedTraversal.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The traversal is defined in the function c08Traversal and uses 8 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID in each direction are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction between particles.
 */
template <class ParticleCell, class Functor>
class C08BasedTraversal : public ColorBasedTraversal<ParticleCell, Functor> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cell block, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit C08BasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor, const double interactionLength,
                             const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : ColorBasedTraversal<ParticleCell, Functor>(dims, functor, interactionLength, cellLength, dataLayout,
                                                   useNewton3) {}

 protected:
  /**
   * The main traversal of the C08Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c08Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class Functor>
template <typename LoopBody>
inline void C08BasedTraversal<ParticleCell, Functor>::c08Traversal(LoopBody &&loopBody) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto end = this->_cellsPerDimension - this->_overlap;
  const auto stride = this->_overlap + 1ul;
  this->colorTraversal(std::forward<LoopBody>(loopBody), end, stride);
}
}  // namespace autopas
