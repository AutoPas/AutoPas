/**
 * @file C18BasedTraversal.h
 * @author nguyen
 * @date 06.09.18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c18 base step.
 *
 * The base step processBaseCell() computes all interactions
 * between the base cell and adjacent cells with greater a ID.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C18BasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C18BasedTraversal(const std::array<unsigned long, 3>& dims, PairwiseFunctor* pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims) {}

  /**
   * C18 traversals are always usable.
   * @return
   */
  bool isApplicable() override { return true; }

 protected:
  /**
   * The main traversal of the C18Traversal. This provides the structure of the loops and its parallelization.
   * @tparam LoopBody
   * @param loopBody the loop body. Normally a lambda function, that takes as as parameters (x,y,z). If you need
   * additional input from outside, please use captures (by reference).
   */
  template <typename LoopBody>
  inline void c18Traversal(LoopBody&& loopBody);
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
template <typename LoopBody>
inline void C18BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::c18Traversal(LoopBody&& loopBody) {
  using std::array;
  const array<unsigned long, 3> stride = {3, 3, 2};
  array<unsigned long, 3> end = {};
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
            loopBody(x, y, z);
          }
        }
      }
    }
  }
}

}  // namespace autopas
