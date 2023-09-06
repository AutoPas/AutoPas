/**
 * @file C04BasedTraversal.h
 * @author C. Menges
 * @date 02.06.2019
 */

#pragma once

#include "ColorBasedTraversal.h"
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
 * @tparam Functor The functor that defines the interaction between particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType, DataLayoutOption::Value dataLayout, bool useNewton3,
          int collapseDepth = 3>
class C04BasedTraversal : public ColorBasedTraversal<ParticleCell, Functor, interactionType, dataLayout, useNewton3, collapseDepth> {
 public:
  /**
   * Constructor of the c04 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length.
   * @param cellLength cell length.
   */
  explicit C04BasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor,
                             const double interactionLength, const std::array<double, 3> &cellLength)
      : ColorBasedTraversal<ParticleCell, Functor, InteractionTypeOption::pairwise, dataLayout, useNewton3, collapseDepth>(
            dims, functor, interactionLength, cellLength) {}

 protected:
  /**
   * The main traversal of the C04Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c04Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType, DataLayoutOption::Value dataLayout, bool useNewton3,
          int collapseDepth>
template <typename LoopBody>
inline void C04BasedTraversal<ParticleCell, Functor, interactionType, dataLayout, useNewton3, collapseDepth>::c04Traversal(
    LoopBody &&loopBody) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto end = this->_cellsPerDimension - this->_overlap;
  auto stride = this->_overlap + 1ul;
  stride[0] = 1;
  this->colorTraversal(std::forward<LoopBody>(loopBody), end, stride);
}
}  // namespace autopas
