/**
 * @file VLCSlicedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellTraversals/SlicedLockBasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCAllCellsNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

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
 * @tparam NeighborList type of the neighbor list
 */

template <class ParticleCell, class PairwiseFunctor, class NeighborList>
class VLCSlicedTraversal : public SlicedLockBasedTraversal<ParticleCell, PairwiseFunctor>,
                           public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param dataLayout
   * @param useNewton3
   * @param typeOfList indicates the type of neighbor list as an enum value, currently only used for getTraversalType
   */
  explicit VLCSlicedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                              double interactionLength, const std::array<double, 3> &cellLength,
                              DataLayoutOption dataLayout, bool useNewton3, ContainerOption::Value typeOfList)
      : SlicedLockBasedTraversal<ParticleCell, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                                dataLayout, useNewton3, false),
        VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList>(typeOfList),
        _functor(pairwiseFunctor) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override {
    switch (this->_typeOfList) {
      case (ContainerOption::verletListsCells):
        return TraversalOption::vlc_sliced;
      case (ContainerOption::pairwiseVerletLists):
        return TraversalOption::vlp_sliced;
      default:
        autopas::utils::ExceptionHandler::exception("Traversal was created with an unsupported neighborlist type: {}",
                                                    this->_typeOfList.to_string());
    }
    // should never be reached.
    return TraversalOption();
  }

  /**
   * VLC Sliced is always applicable to the domain.
   * @return true
   */
  [[nodiscard]] bool isApplicableToDomain() const override { return true; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   * This traversal does not use the CellFunctor, so the function has no effect here
   */
  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, class NeighborList>
inline void VLCSlicedTraversal<ParticleCell, PairwiseFunctor, NeighborList>::traverseParticles() {
  if (this->_dataLayout == DataLayoutOption::soa) {
    this->loadSoA(_functor, *(this->_verletList));
  }

  this->slicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
    auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<PairwiseFunctor>(*(this->_verletList), baseIndex, _functor, this->_dataLayout,
                                                     this->_useNewton3);
  });

  if (this->_dataLayout == DataLayoutOption::soa) {
    this->extractSoA(_functor, *(this->_verletList));
  }
}

}  // namespace autopas
