/**
 * @file VLCCellPairC08Traversal.h
 * @author tirgendetwas
 * @date 31.12.20
 */

#pragma once
#include "autopas/containers/cellPairTraversals/C08BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairC08CellHandler.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * C08 traversal for VLCCellPairNeighborList.
 * The pairwise neighbor list allows access to the relevant pairs of interacting particles for each pair of cells,
 * including the diagonal non-base pair of cells in the standard c08 step.
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>

class VLCCellPairC08Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                                public VLCCellPairTraversalInterface<typename ParticleCell::ParticleType> {
 public:
  /**
   * Constructor of the c08 traversal fot VLCCellPairNeighborList.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   */
  explicit VLCCellPairC08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                   double interactionLength, const std::array<double, 3> &cellLength)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _functor(pairwiseFunctor),
        _cellHandler(dims, pairwiseFunctor, interactionLength, cellLength, this->_overlap) {}

  void traverseParticlePairs() override;

  [[nodiscard]] bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa);
  }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vlp_c08; }

 private:
  PairwiseFunctor *_functor;
  VLCCellPairC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;
  /**
   * Structure of arrays to be used if the data layout is SoA.
   */
  SoA<typename ParticleCell::ParticleType::SoAArraysType> *_soa;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void VLCCellPairC08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  if (dataLayout == DataLayoutOption::soa) {
    _soa = (*(this->_cellPairVerletList)).loadSoA(_functor);
  }

  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    const auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processCellListsC08(*(this->_cellPairVerletList), baseIndex, _functor, dataLayout, _soa,
                                     this->_cellsPerDimension);
  });

  if (dataLayout == DataLayoutOption::soa) {
    (*(this->_cellPairVerletList)).extractSoA(_functor);
    _soa = nullptr;
  }
}

}  // namespace autopas
