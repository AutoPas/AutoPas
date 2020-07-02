/**
 * @file BalancedSlicedTraversalVerlet.h
 *
 * @date 02 May 2020
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellPairTraversals/BalancedSlicedBasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VerletListsCellsTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the balanced sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain along this dimension in one slice (block) per thread. The cut
 * positions are calculated to even out load among the threads. Slices are
 * assigned to the threads in a round robin fashion. Each thread locks the cells
 * on the boundary wall to the previous slice with one lock. This lock is lifted
 * as soon the boundary wall is fully processed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class BalancedSlicedTraversalVerlet
    : public BalancedSlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
      public VerletListsCellsTraversal<typename ParticleCell::ParticleType> {
 public:
  /**
   * Constructor of the balanced sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit BalancedSlicedTraversalVerlet(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                         double interactionLength, const std::array<double, 3> &cellLength)
      : BalancedSlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(
            dims, pairwiseFunctor, interactionLength, cellLength),
        _functor(pairwiseFunctor) {}

  void traverseParticlePairs() override;

  DataLayoutOption getDataLayout() const override { return dataLayout; }

  bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vlc_sliced_balanced; }

  bool isApplicable() const override { return dataLayout == DataLayoutOption::aos; }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void
BalancedSlicedTraversalVerlet<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  this->template slicedTraversal</*allCells*/ true>([&](unsigned long x, unsigned long y, unsigned long z) {
    auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template iterateVerletListsCell<PairwiseFunctor, useNewton3>(*(this->_verletList), baseIndex, _functor);
  });
}

}  // namespace autopas
