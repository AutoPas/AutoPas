/**
 * @file VLCSlicedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellPairTraversals/SlicedBasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into multiple slices along this dimension. Slices are
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
class VLCSlicedTraversal : public SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                           public VLCTraversalInterface<typename ParticleCell::ParticleType> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit VLCSlicedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                              double interactionLength, const std::array<double, 3> &cellLength)
      : SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                    interactionLength, cellLength),
        _functor(pairwiseFunctor) {}

  void traverseParticlePairs() override;

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vlc_sliced; }

  [[nodiscard]] bool isApplicable() const override { return dataLayout == DataLayoutOption::aos; }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void VLCSlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  this->template slicedTraversal</*allCells*/ true>([&](unsigned long x, unsigned long y, unsigned long z) {
    auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template iterateVerletListsCell<PairwiseFunctor, useNewton3>(*(this->_verletList), baseIndex, _functor);
  });
}

}  // namespace autopas
