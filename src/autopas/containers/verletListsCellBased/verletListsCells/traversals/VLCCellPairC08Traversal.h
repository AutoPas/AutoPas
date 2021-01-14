//
// Created by TinaVl on 12/31/2020.
//

#pragma once
#include "autopas/options/DataLayoutOption.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C08BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCCellPairC08CellHandler.h"

namespace autopas
{
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>

class VLCCellPairC08Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                                public VLCCellPairTraversalInterface<typename ParticleCell::ParticleType>
{

 public:
  /**
   * TODO
   */
  explicit VLCCellPairC08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                           double interactionLength, const std::array<double, 3> &cellLength)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _functor(pairwiseFunctor),
        _cellHandler(dims, pairwiseFunctor, interactionLength, cellLength, this->_overlap){}

  void traverseParticlePairs() override;

  [[nodiscard]] bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa);
  }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override {
    return TraversalOption::vlp_c08;
  }

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
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processCellListsC08(*(this->_cellPairVerletList), baseIndex, _functor, dataLayout, _soa, this->_cellsPerDimension);
  });

  if (dataLayout == DataLayoutOption::soa) {
    (*(this->_cellPairVerletList)).extractSoA(_functor);
    _soa = nullptr;
  }
}

}
