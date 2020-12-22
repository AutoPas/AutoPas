/**
 * @file VLCSlicedC02Traversal.h
 *
 * @date 31 May 2020
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellPairTraversals/SlicedC02BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the colored sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into as many slices as possible along this dimension. Unlike the regular
 * sliced traversal, this version uses a 2-coloring to prevent race conditions, instead of
 * locking the starting layers.
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
class VLCSlicedC02Traversal
    : public SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, false>,
      public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
  using VLCTypeOfList = typename VerletListsCellsHelpers<typename ParticleCell::ParticleType>::VLCTypeOfList;

 public:
  /**
   * Constructor of the colored sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit VLCSlicedC02Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                 double interactionLength, const std::array<double, 3> &cellLength)
      : SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, false>(
            dims, pairwiseFunctor, interactionLength, cellLength),
        _functor(pairwiseFunctor) {}

  void traverseParticlePairs() override;

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override {
    switch (typeOfList) {
      case (VLCTypeOfList::vlc):
        return TraversalOption::vlc_sliced_c02;
      case (VLCTypeOfList::vlp):
        return TraversalOption::vlp_sliced_c02;
    }
  }

  [[nodiscard]] bool isApplicable() const override { return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa); }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          class NeighborList,
          typename VerletListsCellsHelpers<typename ParticleCell::ParticleType>::VLCTypeOfList::Value typeOfList>
inline void VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, NeighborList,
                                  typeOfList>::traverseParticlePairs() {
  if(dataLayout == DataLayoutOption::soa)
  {
    this->setupLoadSoA(_functor);
  }

  this->cSlicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<PairwiseFunctor, useNewton3>(*(this->_verletList), baseIndex, _functor, dataLayout);
  });

  if(dataLayout == DataLayoutOption::soa)
  {
    this->setupExtractSoA(_functor);
  }
}

}  // namespace autopas
