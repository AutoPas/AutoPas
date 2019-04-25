/**
 * @file C01BasedTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <autopas/utils/WrapOpenMP.h>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c01 base step.
 *
 * The traversal is defined in the function c01Traversal and uses 1 color. Interactions between two cells are allowed
 * only if particles of the first cell are modified. This means that newton3 optimizations are NOT allowed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C01BasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit C01BasedTraversal(const std::array<unsigned long, 3>& dims, PairwiseFunctor* pairwiseFunctor,
                             double cutoff = 1.0, const std::array<double, 3>& cellLength = {1.0, 1.0, 1.0})
      : CellPairTraversal<ParticleCell>(dims), _cutoff(cutoff), _cellLength(cellLength) {
    for (unsigned int d = 0; d < 3; d++) {
      _overlap[d] = std::ceil(_cutoff / _cellLength[d]);
    }
  }

  /**
   * C01 traversals are only usable if useNewton3 is disabled.
   *
   * This is because the cell functor in the c01 traversal is hardcoded to not allow newton 3 even if only one thread is
   * used.
   *
   * @return
   */
  bool isApplicable() override { return not useNewton3; }

 protected:
  /**
   * The main traversal of the C01Traversal.
   * This provides the structure of the loops and its parallelization.
   * @tparam LoopBody
   * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
   * (x,y,z). If you need additional input from outside, please use captures (by reference).
   */
  template <typename LoopBody>
  inline void c01Traversal(LoopBody&& loopBody);

  /**
   * cutoff radius.
   */
  double _cutoff;

  /**
   * cell length in CellBlock3D.
   */
  std::array<double, 3> _cellLength;

  /**
   * overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  std::array<unsigned long, 3> _overlap;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
template <typename LoopBody>
inline void C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::c01Traversal(LoopBody&& loopBody) {
  const unsigned long end_x = this->_cellsPerDimension[0] - _overlap[0];
  const unsigned long end_y = this->_cellsPerDimension[1] - _overlap[1];
  const unsigned long end_z = this->_cellsPerDimension[2] - _overlap[2];

#if defined(AUTOPAS_OPENMP)
  // @todo: find optimal chunksize
#pragma omp parallel for schedule(dynamic) collapse(3)
#endif
  for (unsigned long z = _overlap[0]; z < end_z; ++z) {
    for (unsigned long y = _overlap[1]; y < end_y; ++y) {
      for (unsigned long x = _overlap[2]; x < end_x; ++x) {
        loopBody(x, y, z);
      }
    }
  }
}

}  // namespace autopas
