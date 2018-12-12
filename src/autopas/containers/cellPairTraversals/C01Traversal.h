/**
 * @file C01Traversal.h
 * @author nguyen
 * @date 16.09.2018
 */

#pragma once

#include "autopas/containers/VerletListsCellsHelpers.h"
#include "autopas/containers/cellPairTraversals/C01BasedTraversal.h"
#include "autopas/containers/cellPairTraversals/VerletListsCellsTraversal.h"
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
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA>
class C01Traversal : public C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>,
                     public VerletListsCellsTraversal<PairwiseFunctor, false> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>(dims, pairwiseFunctor),
        VerletListsCellsTraversal<PairwiseFunctor, false>(pairwiseFunctor) {}
  // documentation in base class
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  /**
   * Traverse verlet lists of all cells.
   * This function needs to be implemented by derived classes.
   * @param verlet verlet lists for each cell
   */
  template <class Particle>
  void traverseCellVerlet(std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>> &verlet);
  TraversalOptions getTraversalType() override;
  bool isApplicable() override;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
inline bool C01Traversal<ParticleCell, PairwiseFunctor, useSoA>::isApplicable() {
  return true;
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, useSoA>::traverseCellPairs(std::vector<ParticleCell> &cells) {
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

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
template <class Particle>
inline void C01Traversal<ParticleCell, PairwiseFunctor, useSoA>::traverseCellVerlet(
    std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>> &verlet) {
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

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
TraversalOptions C01Traversal<ParticleCell, PairwiseFunctor, useSoA>::getTraversalType() {
  return TraversalOptions::c01;
};

}  // namespace autopas
