/**
 * @file VLCC08Traversal.h
 * @date 16 May 2024
 * @author F.Gratl
 */

#pragma once

#include "autopas/containers/cellTraversals/C08BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"

namespace autopas {

/**
 * This class provides the vlc_c08 traversal.
 *
 * The traversal uses the c08 base step performed on every single cell. Since
 * these steps overlap a domain coloring with eight colors is applied.
 *
 * For each cell all neighbor lists are processed, so depending on whether lists
 * were built with newton3 the base step is c01 or c08
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of two particles.
 * @tparam NeighborList type of the neighbor list
 */
template <class ParticleCell, class Functor, class NeighborList>
class VLCC08Traversal : public C08BasedTraversal<ParticleCell, Functor>,
                        public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
 public:
  /**
   * Constructor of the vlc_c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLCC08Traversal(const std::array<unsigned long, 3> &dims, Functor *pairwiseFunctor, double interactionLength,
                           const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : C08BasedTraversal<ParticleCell, Functor>(dims, pairwiseFunctor, interactionLength, cellLength, dataLayout,
                                                 useNewton3),
        VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList>(ContainerOption::verletListsCells),
        _functor(pairwiseFunctor) {}

  void traverseParticles() override;

  void traverseParticlePairs();
  void traverseParticleTriplets();

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vlc_c08; }

  [[nodiscard]] bool isApplicable() const override {
    // This traversal is only safe to use for CSF>=1
    const double minCellLength = *std::min_element(this->_cellLength.cbegin(), this->_cellLength.cend());
    const bool maxOneCellInCutoff = minCellLength >= this->_interactionLength;

    return maxOneCellInCutoff;
  }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   * This traversal does not use the CellFunctor, so the function has no effect here
   */
  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  Functor *_functor;
};

template <class ParticleCell, class Functor, class NeighborList>
inline void VLCC08Traversal<ParticleCell, Functor, NeighborList>::traverseParticles() {
  if constexpr (utils::isPairwiseFunctor<Functor>()) {
    traverseParticlePairs();
  } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
    traverseParticleTriplets();
  } else {
    utils::ExceptionHandler::exception(
        "VLListIterationTraversal::traverseParticles(): Functor {} is not of type PairwiseFunctor or TriwiseFunctor.",
        _functor->getName());
  }
}

template <class ParticleCell, class Functor, class NeighborList>
inline void VLCC08Traversal<ParticleCell, Functor, NeighborList>::traverseParticlePairs() {
  if (this->_dataLayout == DataLayoutOption::soa) {
    this->loadSoA(_functor, *(this->_verletList));
  }

  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    const auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<Functor>(*(this->_verletList), baseIndex, _functor, this->_dataLayout,
                                             this->_useNewton3);
  });

  if (this->_dataLayout == DataLayoutOption::soa) {
    this->extractSoA(_functor, *(this->_verletList));
  }
}

template <class ParticleCell, class Functor, class NeighborList>
inline void VLCC08Traversal<ParticleCell, Functor, NeighborList>::traverseParticleTriplets() {
  if (this->_dataLayout == DataLayoutOption::soa) {
    this->loadSoA(_functor, *(this->_verletList));
  }

  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    const auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<Functor>(*(this->_verletList), baseIndex, _functor, this->_dataLayout,
                                             this->_useNewton3);
  });

  if (this->_dataLayout == DataLayoutOption::soa) {
    this->extractSoA(_functor, *(this->_verletList));
  }
}

}  // namespace autopas