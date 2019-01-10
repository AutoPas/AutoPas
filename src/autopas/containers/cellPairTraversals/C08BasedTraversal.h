/**
 * @file C08BasedTraversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The base step processBaseCell() computes one set of pairwise interactions
 * between two cells for each spatial direction based on the baseIndex.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C08BasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C08BasedTraversal(const std::array<unsigned long, 3>& dims, PairwiseFunctor* pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims) {}

  /**
   * C08 traversals are always usable.
   * @return
   */
  bool isApplicable() override { return true; }

 protected:
  /**
   * The main traversal of the C08Traversal. This provides the structure of the loops and its parallelization.
   * @tparam LoopBody
   * @param loopBody the loop body. Normally a lambda function, that takes as as parameters (x,y,z). If you need
   * additional input from outside, please use captures (by reference).
   */
  template <typename LoopBody>
  inline void c08Traversal(LoopBody&& loopBody);
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
template <typename LoopBody>
inline void C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::c08Traversal(LoopBody&& loopBody) {
  using std::array;
  const array<unsigned long, 3> stride = {2, 2, 2};
  array<unsigned long, 3> end = {};
  for (int d = 0; d < 3; ++d) {
    end[d] = this->_cellsPerDimension[d] - 1;
  }

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
    for (unsigned long col = 0; col < 8; ++col) {
      std::array<unsigned long, 3> start = utils::ThreeDimensionalMapping::oneToThreeD(col, stride);

      // intel compiler demands following:
      const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
      const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
      const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3)
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
