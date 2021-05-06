/**
 * @file VLCC18Traversal.h
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include "autopas/containers/cellPairTraversals/C18BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/WrapOpenMP.h"

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
 * @tparam dataLayout
 * @tparam useNewton3
 * @tparam NeighborList type of the neighbor list
 * @tparam typeOfList indicates the type of neighbor list as an enum value, currently only used for getTraversalType
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          class NeighborList, ContainerOption::Value typeOfList>
class VLCC18Traversal : public C18BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                        public VLCTraversalInterface<typename ParticleCell::ParticleType, NeighborList> {
 public:
  /**
   * Constructor of the lc_c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param interactionLength cutoff + skin
   * @param cellLength length of the underlying cells
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit VLCC18Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                           const double interactionLength, const std::array<double, 3> &cellLength)
      : C18BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _functor(pairwiseFunctor) {}

  void traverseParticlePairs() override;

  [[nodiscard]] TraversalOption getTraversalType() const override {
    switch (typeOfList) {
      case (ContainerOption::verletListsCells):
        return TraversalOption::vlc_c18;
      case (ContainerOption::pairwiseVerletLists):
        return TraversalOption::vlp_c18;
    }
    // should never be reached.
    return TraversalOption();
  }

  [[nodiscard]] bool isApplicable() const override { return dataLayout == DataLayoutOption::aos; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          class NeighborList, ContainerOption::Value typeOfList>
inline void VLCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, NeighborList,
                            typeOfList>::traverseParticlePairs() {
  this->template c18Traversal</*allCells*/ true>([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->template processCellLists<PairwiseFunctor, useNewton3>(*(this->_verletList), baseIndex, _functor);
  });
}

}  // namespace autopas
