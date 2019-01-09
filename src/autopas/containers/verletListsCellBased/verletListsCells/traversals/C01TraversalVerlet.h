/**
 * @file C01TraversalVerlet.h
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include "autopas/containers/cellPairTraversals/C01BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VerletListsCellsTraversal.h"
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
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C01TraversalVerlet
    : public C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>,
      public VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C01TraversalVerlet(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>(dims, pairwiseFunctor),
        VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3>(pairwiseFunctor) {}
  // documentation in base class
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  /**
   * @copydoc VerletListsCellsTraversal::traverseCellVerlet
   */
  void traverseCellVerlet(typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                                             useNewton3>::verlet_storage_type &verlet) override;
  TraversalOptions getTraversalType() override {return TraversalOptions::c01;}
  bool isApplicable() override {return not useNewton3;}
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C01TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellVerlet(
    typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                       useNewton3>::verlet_storage_type &verlet) {
  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    this->iterateVerletListsCell(verlet, baseIndex);
  });
}

}  // namespace autopas
