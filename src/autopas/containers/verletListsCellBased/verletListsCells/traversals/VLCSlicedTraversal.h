/**
 * @file VLCSlicedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellPairTraversals/SlicedLockBasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the (locked) sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into multiple slices along this dimension. Slices are
 * assigned to the threads in a round robin fashion. Each thread locks the cells
 * on the boundary wall to the previous slice with one lock. This lock is lifted
 * as soon the boundary wall is fully processed.
 *
 * For each cell all neighbor lists are processed, so depending on whether lists
 * were built with newton3 the base step is c01 or c18
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 * @tparam NeighborList type of the neighbor list
 * @tparam typeOfList indicates the type of neighbor list as an enum value, currently only used for getTraversalType
 */

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          class NeighborList,
          typename VerletListsCellsHelpers<typename ParticleCell::ParticleType>::VLCTypeOfList::Value typeOfList>
class VLCSlicedTraversal
    : public SlicedLockBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, false>,
      public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
  using VLCTypeOfList = typename VerletListsCellsHelpers<typename ParticleCell::ParticleType>::VLCTypeOfList;

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
      : SlicedLockBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, false>(
            dims, pairwiseFunctor, interactionLength, cellLength),
        _functor(pairwiseFunctor) {}

  void traverseParticlePairs() override;

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override {
    switch (typeOfList) {
      case (VLCTypeOfList::vlc):
        return TraversalOption::vlc_sliced;
      case (VLCTypeOfList::vlp):
        return TraversalOption::vlp_sliced;
    }
  }

  [[nodiscard]] bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa);
  }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          class NeighborList,
          typename VerletListsCellsHelpers<typename ParticleCell::ParticleType>::VLCTypeOfList::Value typeOfList>
inline void VLCSlicedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, NeighborList,
                               typeOfList>::traverseParticlePairs() {
  if (dataLayout == DataLayoutOption::soa) {
    this->setupLoadSoA(_functor, *(this->_verletList));
  }

  this->slicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<PairwiseFunctor, useNewton3>(*(this->_verletList), baseIndex, _functor, dataLayout);
  });

  if (dataLayout == DataLayoutOption::soa) {
    this->setupExtractSoA(_functor, *(this->_verletList));
  }
}

}  // namespace autopas
