/**
 * @file C01BasedTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include "ColorBasedTraversal.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c01 base step.
 *
 * The traversal is defined in the function c01Traversal and uses 1 color. Interactions between cells are allowed
 * only if particles of only one cell are modified. This means that newton3 optimizations are NOT allowed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of particles.
 * @tparam dataLayout indicates usage of SoA
 */
template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType,
          int collapseDepth = 3>
class C01BasedTraversal
    : public ColorBasedTraversal<ParticleCell, Functor, interactionType, collapseDepth> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit C01BasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor, double interactionLength,
                             const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : ColorBasedTraversal<ParticleCell, Functor, interactionType, collapseDepth>(
            dims, functor, interactionLength, cellLength, dataLayout, useNewton3) {}

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

template <class ParticleCell, class Functor, InteractionTypeOption::Value interactionType, int collapseDepth>
template <typename LoopBody>
inline void C01BasedTraversal<ParticleCell, Functor, interactionType,
                              collapseDepth>::c01Traversal(LoopBody &&loopBody) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto offset = this->_overlap;
  const auto end = this->_cellsPerDimension - this->_overlap;
  this->colorTraversal(std::forward<LoopBody>(loopBody), end, {1ul, 1ul, 1ul}, offset);
}
}  // namespace autopas
