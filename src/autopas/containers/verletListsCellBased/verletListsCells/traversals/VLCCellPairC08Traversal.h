/**
 * @file VLCCellPairC08Traversal.h
 * @author tirgendetwas
 * @date 31.12.20
 */

#pragma once
#include "autopas/containers/cellTraversals/C08BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairC08CellHandler.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * C08 traversal for VLCCellPairNeighborList.
 * The pairwise neighbor list allows access to the relevant pairs of interacting particles for each pair of cells,
 * including the diagonal non-base pair of cells in the standard c08 step.
 * @tparam ParticleCell_T the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor>
class VLCCellPairC08Traversal : public C08BasedTraversal<ParticleCell_T, PairwiseFunctor>,
                                public VLCCellPairTraversalInterface<typename ParticleCell_T::ParticleType> {
 public:
  /**
   * Constructor of the c08 traversal for VLCCellPairNeighborList.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLCCellPairC08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor &pairwiseFunctor,
                                   double interactionLength, const std::array<double, 3> &cellLength,
                                   DataLayoutOption dataLayout, bool useNewton3)
      : C08BasedTraversal<ParticleCell_T, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                           dataLayout, useNewton3),
        _functor(pairwiseFunctor),
        _cellHandler(dims, interactionLength, cellLength) {}

  void traverseParticles() override;

  [[nodiscard]] bool isApplicable() const override {
    return (this->_dataLayout == DataLayoutOption::aos or this->_dataLayout == DataLayoutOption::soa);
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vlp_c08; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   * This traversal does not use the CellFunctor, so the function has no effect here
   */
  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  PairwiseFunctor &_functor;
  VLCCellPairC08CellHandler<ParticleCell_T, PairwiseFunctor> _cellHandler;
  /**
   * Structure of arrays to be used if the data layout is SoA.
   */
  SoA<typename ParticleCell_T::ParticleType::SoAArraysType> *_soa;
};

template <class ParticleCell_T, class PairwiseFunctor>
inline void VLCCellPairC08Traversal<ParticleCell_T, PairwiseFunctor>::traverseParticles() {
  if (this->_dataLayout == DataLayoutOption::soa) {
    _soa = this->_cellPairVerletList->loadSoA(_functor);
  }

  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    const auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processCellListsC08(*(this->_cellPairVerletList), baseIndex, _functor, this->_dataLayout, _soa,
                                     this->_useNewton3);
  });

  if (this->_dataLayout == DataLayoutOption::soa) {
    this->_cellPairVerletList->extractSoA(_functor);
    _soa = nullptr;
  }
}

}  // namespace autopas
