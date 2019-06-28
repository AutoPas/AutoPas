/**
 * @file C18TraversalVerlet.h
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include "autopas/containers/cellPairTraversals/C18BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VerletListsCellsTraversal.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c18 traversal.
 *
 * The traversal uses the c18 base step performed on every single cell. Since
 * these steps overlap a domain coloring with eighteen colors is applied.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C18TraversalVerlet
    : public C18BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>,
      public VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3> {
 public:
  /**
   * Constructor of the c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C18TraversalVerlet(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : C18BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor),
        VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3>(pairwiseFunctor) {}

  /**
   * @copydoc VerletListsCellsTraversal::traverseCellVerlet
   */
  void traverseCellVerlet(typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                                             useNewton3>::verlet_storage_type &verlet) override;

  TraversalOption getTraversalType() const override { return TraversalOption::c18Verlet; };

  bool isApplicable() const override { return DataLayout == DataLayoutOption::aos; }
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C18TraversalVerlet<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellVerlet(
    typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                       useNewton3>::verlet_storage_type &verlet) {
  this->c18Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->iterateVerletListsCell(verlet, baseIndex);
  });
}

}  // namespace autopas
