/**
 * @file VLCSlicedC02Traversal.h
 *
 * @date 31 May 2020
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellTraversals/SlicedC02BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

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
 * @tparam NeighborList type of the neighbor list
 */
template <class ParticleCell, class PairwiseFunctor, class NeighborList>
class VLCSlicedC02Traversal : public SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor>,
                              public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
 public:
  /**
   * Constructor of the colored sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param dataLayout
   * @param useNewton3
   * @param typeOfList indicates the type of neighbor list as an enum value, currently only used for getTraversalType
   */
  explicit VLCSlicedC02Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                 double interactionLength, const std::array<double, 3> &cellLength,
                                 DataLayoutOption dataLayout, bool useNewton3, ContainerOption::Value typeOfList)
      : SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                               dataLayout, useNewton3, false),
        VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList>(typeOfList),
        _functor(pairwiseFunctor) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override {
    switch (this->_typeOfList) {
      case (ContainerOption::verletListsCells):
        return TraversalOption::vlc_sliced_c02;
      case (ContainerOption::pairwiseVerletLists):
        return TraversalOption::vlp_sliced_c02;
      default:
        autopas::utils::ExceptionHandler::exception("Traversal was created with an unsupported neighborlist type: {}",
                                                    this->_typeOfList.to_string());
    }
    // should never be reached.
    return TraversalOption();
  }

  [[nodiscard]] bool isApplicable() const override {
    return (this->_dataLayout == DataLayoutOption::aos or this->_dataLayout == DataLayoutOption::soa);
  }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   * This traversal does not use the CellFunctor, so the function has no effect here
   */
  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, class NeighborList>
inline void VLCSlicedC02Traversal<ParticleCell, PairwiseFunctor, NeighborList>::traverseParticles() {
  if (this->_dataLayout == DataLayoutOption::soa) {
    this->loadSoA(_functor, *(this->_verletList));
  }

  this->cSlicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<PairwiseFunctor>(*(this->_verletList), baseIndex, _functor, this->_dataLayout,
                                                     this->_useNewton3);
  });

  if (this->_dataLayout == DataLayoutOption::soa) {
    this->extractSoA(_functor, *(this->_verletList));
  }
}

}  // namespace autopas
