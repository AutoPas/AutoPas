/**
 * @file C08BasedTraversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/BlackBoxOptions.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The traversal is defined in the function c08Traversal and uses 8 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID in each direction are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 * @tparam blackBoxTraversalOption The blackbox traversal option to be used.
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3,
          BlackBoxTraversalOption blackBoxTraversalOption = normal>
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
   * The main traversal of the C08Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c08Traversal(LoopBody&& loopBody);

  /**
   * The main traversal of the C08Traversal for the outer blackbox case.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void c08TraversalOuter(LoopBody&& loopBody);
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3,
          BlackBoxTraversalOption blackBoxTraversalOption>
template <typename LoopBody>
inline void C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3, blackBoxTraversalOption>::c08Traversal(
    LoopBody&& loopBody) {
  if (blackBoxTraversalOption == outer) {
    // tooSmall indicates whether there is any inner region.
    // If there is no inner region, just do a normal traversal, as the traversals are supposedly the same.
    bool tooSmall = false;
    for (unsigned short d = 0; d < 3; ++d) {
      tooSmall |= this->_cellsPerDimension[d] <= 6;
    }
    if (not tooSmall) {
      return c08TraversalOuter(loopBody);
    }  // else: do normal traversal!
  }
  using std::array;
  const array<unsigned long, 3> stride = {2, 2, 2};
  const array<unsigned long, 3> blackBoxOffset = {2, 2, 2};
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

      if (blackBoxTraversalOption == inner) {
        start = ArrayMath::add(start, blackBoxOffset);
        end = ArrayMath::sub(end, blackBoxOffset);
      }

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

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3,
          BlackBoxTraversalOption blackBoxTraversalOption>
template <typename LoopBody>
inline void C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3,
                              blackBoxTraversalOption>::c08TraversalOuter(LoopBody&& loopBody) {
  using std::array;
  const array<unsigned long, 3> stride = {2, 2, 2};
  array<unsigned long, 3> end = {};
  for (int d = 0; d < 3; ++d) {
    end[d] = this->_cellsPerDimension[d];
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

      // z:
      {
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3) nowait
#endif
        // lower z
        for (unsigned long z = start_z; z < 3; z += stride_z) {
          for (unsigned long y = start_y; y < end_y; y += stride_y) {
            for (unsigned long x = start_x; x < end_x; x += stride_x) {
              loopBody(x, y, z);
            }
          }
        }
      }

      unsigned long end_z_inner;
      {
        // upper z: always only one value for z:
        unsigned long z;
        if (start_z == 1) {
          // round to first uneven number smaller than end_z
          z = (end_z & ~1ul) - 1ul;
        } else {  // start_z == 0
          // round to first even number smaller than end_z
          z = (end_z - 1ul) & ~1ul;
        }
        end_z_inner = z;
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(2) nowait
#endif
        for (unsigned long y = start_y; y < end_y; y += stride_y) {
          for (unsigned long x = start_x; x < end_x; x += stride_x) {
            loopBody(x, y, z);
          }
        }
      }

      const unsigned long start_z_inner = 4 - start_z;  // (start_z == 0 ? 4 : 3);
      // y:
      {
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3) nowait
#endif
        // lower y
        for (unsigned long z = start_z_inner; z < end_z_inner; z += stride_z) {
          for (unsigned long y = start_y; y < 3; y += stride_y) {
            for (unsigned long x = start_x; x < end_x; x += stride_x) {
              loopBody(x, y, z);
            }
          }
        }
      }

      unsigned long end_y_inner;
      {
        // upper y: always only one value for y:
        unsigned long y;
        if (start_y % 2 == 1) {
          // round to first uneven number smaller than end_z
          y = (end_y & ~1ul) - 1ul;
        } else {  // start_z == 0
          // round to first even number smaller than end_z
          y = (end_y - 1ul) & ~1ul;
        }
        end_y_inner = y;
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(2) nowait
#endif
        for (unsigned long z = start_z_inner; z < end_z_inner; z += stride_z) {
          for (unsigned long x = start_x; x < end_x; x += stride_x) {
            loopBody(x, y, z);
          }
        }
      }

      const unsigned long start_y_inner = 4 - start_y;  // (start_z == 0 ? 4 : 3);

      // x:
      {
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3) nowait
#endif
        // lower y
        for (unsigned long z = start_z_inner; z < end_z_inner; z += stride_z) {
          for (unsigned long y = start_y_inner; y < end_y_inner; y += stride_y) {
            for (unsigned long x = start_x; x < 3; x += stride_x) {
              loopBody(x, y, z);
            }
          }
        }
      }

      {
        // upper y: always only one value for y:
        unsigned long x;
        if (start_x % 2 == 1) {
          // round to first uneven number smaller than end_z
          x = (end_x & ~1ul) - 1ul;
        } else {  // start_z == 0
          // round to first even number smaller than end_z
          x = (end_x - 1ul) & ~1ul;
        }
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(2)
#endif
        for (unsigned long z = start_z_inner; z < end_z_inner; z += stride_z) {
          for (unsigned long y = start_y_inner; y < end_y_inner; y += stride_y) {
            loopBody(x, y, z);
          }
        }
      }
    }
  }
}

}  // namespace autopas
