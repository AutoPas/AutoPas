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
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C18TraversalVerlet
    : public C18BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>,
      public VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3> {
 public:
  /**
   * Constructor of the c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C18TraversalVerlet(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : C18BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>(dims, pairwiseFunctor),
        VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor, useNewton3>(pairwiseFunctor) {}

  /**
   * @copydoc VerletListsCellsTraversal::traverseCellVerlet
   */
  void traverseCellVerlet(typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                                             useNewton3>::verlet_storage_type &verlet) override;
  TraversalOptions getTraversalType() override;
  bool isApplicable() override;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline bool C18TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::isApplicable() {
  return true;
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C18TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellVerlet(
    typename VerletListsCellsTraversal<typename ParticleCell::ParticleType, PairwiseFunctor,
                                       useNewton3>::verlet_storage_type &verlet) {
  using std::array;
  const array<unsigned long, 3> stride = {3, 3, 2};
  array<unsigned long, 3> end;
  end[0] = this->_cellsPerDimension[0];
  end[1] = this->_cellsPerDimension[1];
  end[2] = this->_cellsPerDimension[2] - 1;

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
    for (unsigned long col = 0; col < 18; ++col) {
      std::array<unsigned long, 3> start = utils::ThreeDimensionalMapping::oneToThreeD(col, stride);

      // intel compiler demands following:
      const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
      const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
      const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

#if defined(AUTOPAS_OPENMP)
      // @todo: find optimal chunksize
#pragma omp for schedule(dynamic) collapse(3)
#endif
      for (unsigned long z = start_z; z < end_z; z += stride_z) {
        for (unsigned long y = start_y; y < end_y; y += stride_y) {
          for (unsigned long x = start_x; x < end_x; x += stride_x) {
            unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
            this->iterateVerletListsCell(verlet, baseIndex);
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
TraversalOptions C18TraversalVerlet<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::getTraversalType() {
  return TraversalOptions::c18;
};

}  // namespace autopas
