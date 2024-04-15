/**
 * @file C08BasedTraversalListBuilding3B.h
 * @author Alexander Haberl
 * @date 2024/04/15
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
 * @tparam PairwiseFunctor The functor that defines the interaction between particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType,
          DataLayoutOption::Value dataLayout, bool useNewton3>
class C08BasedNeighborListBuilding3B
    : public ColorBasedTraversal<ParticleCell, Functor, interactionType, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cell block, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit C08BasedNeighborListBuilding3B(const std::array<unsigned long, 3> &dims, Functor *functor,
                                          const double interactionLength, const std::array<double, 3> &cellLength)
      : ColorBasedTraversal<ParticleCell, Functor, interactionType, dataLayout, useNewton3>(
            dims, functor, interactionLength, cellLength) {}

 protected:
  /**
   * The main traversal of the C08Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c08Traversal(LoopBody &&loopBody);
};

template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType,
          DataLayoutOption::Value dataLayout, bool useNewton3>
template <typename LoopBody>
inline void C08BasedNeighborListBuilding3B<ParticleCell, Functor, interactionType, dataLayout,
                                           useNewton3>::c08Traversal(LoopBody &&loopBody) {
  using namespace autopas::utils::ArrayMath::literals;

  // last cells also have to be traversed to have correct Neighborlist building between Halo Cells
  const auto end = this->_cellsPerDimension;
  const auto stride = this->_overlap + 1ul;
  this->colorTraversal(std::forward<LoopBody>(loopBody), end, stride);
}
}  // namespace autopas
