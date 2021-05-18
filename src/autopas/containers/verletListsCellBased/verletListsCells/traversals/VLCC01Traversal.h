/**
 * @file VLCC01Traversal.h
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include "autopas/containers/cellPairTraversals/C01BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c01 traversal.
 *
 * The traversal uses the c01 base step performed on every single cell.
 * newton3 cannot be applied!
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 * @tparam NeighborList type of the neighbor list
 * @tparam typeOfList indicates the type of neighbor list as an enum value, currently only used for getTraversalType
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          class NeighborList, ContainerOption::Value typeOfList>
class VLCC01Traversal : public C01BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                        public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit VLCC01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                           double interactionLength, const std::array<double, 3> &cellLength)
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _functor(pairwiseFunctor) {}

  void traverseParticlePairs() override;

  [[nodiscard]] TraversalOption getTraversalType() const override {
    switch (typeOfList) {
      case (ContainerOption::verletListsCells):
        return TraversalOption::vlc_c01;
      case (ContainerOption::pairwiseVerletLists):
        return TraversalOption::vlp_c01;
    }
    // should never be reached.
    return TraversalOption();
  }

  [[nodiscard]] bool isApplicable() const override {
    return (not useNewton3) and (dataLayout == DataLayoutOption::aos or dataLayout == DataLayoutOption::soa);
  }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          class NeighborList, ContainerOption::Value typeOfList>
inline void VLCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, NeighborList,
                            typeOfList>::traverseParticlePairs() {
  if (dataLayout == DataLayoutOption::soa) {
    this->loadSoA(_functor, *(this->_verletList));
  }

  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<PairwiseFunctor, useNewton3>(*(this->_verletList), baseIndex, _functor, dataLayout);
  });

  if (dataLayout == DataLayoutOption::soa) {
    this->extractSoA(_functor, *(this->_verletList));
  }
}

}  // namespace autopas
