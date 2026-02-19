/**
 * @file VLCC18Traversal.h
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include "autopas/containers/cellTraversals/C18BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"

namespace autopas {

/**
 * This class provides the lc_c18 traversal.
 *
 * The traversal uses the c18 base step performed on every single cell. Since
 * these steps overlap a domain coloring with eighteen colors is applied.
 *
 * For each cell all neighbor lists are processed, so depending on whether lists
 * were built with newton3 the base step is c01 or c18
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam NeighborList type of the neighbor list
 */
template <class ParticleCell, class PairwiseFunctor, class NeighborList>
class VLCC18Traversal : public C18BasedTraversal<ParticleCell, PairwiseFunctor>,
                        public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
 public:
  /**
   * Constructor of the lc_c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param dataLayout
   * @param useNewton3
   * @param typeOfList indicates the type of neighbor list as an enum value, currently only used for getTraversalType
   */
  explicit VLCC18Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                           const double interactionLength, const std::array<double, 3> &cellLength,
                           DataLayoutOption dataLayout, bool useNewton3, ContainerOption::Value typeOfList)
      : C18BasedTraversal<ParticleCell, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                         dataLayout, useNewton3),
        VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList>(typeOfList),
        _functor(pairwiseFunctor) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override {
    switch (this->_typeOfList) {
      case (ContainerOption::verletListsCells):
        return TraversalOption::vlc_c18;
      case (ContainerOption::pairwiseVerletLists):
        return TraversalOption::vlp_c18;
      default:
        utils::ExceptionHandler::exception("Traversal was created with an unsupported neighborlist type: {}",
                                           this->_typeOfList.to_string());
    }
    // should never be reached.
    return TraversalOption();
  }

  /**
   * VLC C18 is always applicable to the domain.
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
void VLCC18Traversal<ParticleCell, PairwiseFunctor, NeighborList>::traverseParticles() {
  if (this->_dataLayout == DataLayoutOption::soa) {
    this->loadSoA(_functor, *(this->_verletList));
  }

  this->template c18Traversal</*allCells*/ true>([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<PairwiseFunctor>(*(this->_verletList), baseIndex, _functor, this->_dataLayout,
                                                     this->_useNewton3);
  });

  if (this->_dataLayout == DataLayoutOption::soa) {
    this->extractSoA(_functor, *(this->_verletList));
  }
}

}  // namespace autopas