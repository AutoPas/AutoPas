/**
 * @file C01Traversal.h
 * @author nguyen
 * @date 16.09.2018
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
class C01Traversal
    : public C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>,
      public VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>(dims, pairwiseFunctor),
        VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3>(pairwiseFunctor) {}
  // documentation in base class
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  /**
   * @copydoc VerletListsCellsTraversal::traverseCellVerlet
   */
  void traverseCellVerlet(typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                                             useNewton3>::verlet_storage_type &verlet) override;
  TraversalOptions getTraversalType() override;
  bool isApplicable() override;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline bool C01Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::isApplicable() {
  return not useNewton3;
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C01Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  const unsigned long end_x = this->_cellsPerDimension[0] - 1;
  const unsigned long end_y = this->_cellsPerDimension[1] - 1;
  const unsigned long end_z = this->_cellsPerDimension[2] - 1;

#if defined(AUTOPAS_OPENMP)
  // @todo: find optimal chunksize
#pragma omp parallel for schedule(dynamic) collapse(3)
#endif
  for (unsigned long z = 1; z < end_z; ++z) {
    for (unsigned long y = 1; y < end_y; ++y) {
      for (unsigned long x = 1; x < end_x; ++x) {
        this->processBaseCell(cells, x, y, z);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C01Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellVerlet(
    typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                       useNewton3>::verlet_storage_type &verlet) {
  const unsigned long end_x = this->_cellsPerDimension[0] - 1;
  const unsigned long end_y = this->_cellsPerDimension[1] - 1;
  const unsigned long end_z = this->_cellsPerDimension[2] - 1;

#if defined(AUTOPAS_OPENMP)
  // @todo: find optimal chunksize
#pragma omp parallel for schedule(dynamic) collapse(3)
#endif
  for (unsigned long z = 1; z < end_z; ++z) {
    for (unsigned long y = 1; y < end_y; ++y) {
      for (unsigned long x = 1; x < end_x; ++x) {
        unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
        this->iterateVerletListsCell(verlet, baseIndex);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
TraversalOptions C01Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::getTraversalType() {
  return TraversalOptions::c01;
};

}  // namespace autopas
