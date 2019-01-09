/**
 * @file C01BasedTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c01 base step.
 *
 * The base step processBaseCell() computes all interactions
 * between the base cell and adjacent cells.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA>
class C01BasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C01BasedTraversal(const std::array<unsigned long, 3>& dims, PairwiseFunctor* pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims),
        _cellFunctor(
            CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, false, false>(
                pairwiseFunctor)) {
    computeOffsets();
  }

 protected:
  template <typename LoopBody>
  inline void c01Traversal(LoopBody&& loopBody);

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, false, false> _cellFunctor;

  /**
   * Pairs for processBaseCell().
   */
  std::vector<int> _cellOffsets;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
inline void C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>::computeOffsets() {
  for (int z = -1; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        int offset = (z * this->_cellsPerDimension[1] + y) * this->_cellsPerDimension[0] + x;
        _cellOffsets.push_back(offset);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
template <typename LoopBody>
inline void C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>::c01Traversal(LoopBody&& loopBody) {
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
        loopBody(x, y, z);
      }
    }
  }
}

}  // namespace autopas
