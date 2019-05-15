/**
 * @file CBasedTraversal.h
 * @author C. Menges
 * @date 26.04.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using base steps based on cell coloring.
 *
 * @tparam ParticleCell the type of cells
 */
template <class ParticleCell>
class CBasedTraversal : public CellPairTraversal<ParticleCell> {
  using ParticleFloatType = typename ParticleCell::ParticleType::ParticleFloatingPointType;

 protected:
  /**
   * Constructor of the CBasedTraversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit CBasedTraversal(const std::array<unsigned long, 3>& dims, const ParticleFloatType cutoff,
                           const std::array<ParticleFloatType, 3>& cellLength)
      : CellPairTraversal<ParticleCell>(dims), _cutoff(cutoff), _cellLength(cellLength) {
    for (unsigned int d = 0; d < 3; d++) {
      _overlap[d] = std::ceil(_cutoff / _cellLength[d]);
    }
  }

  /**
   * Destructor of CBasedTraversal.
   */
  ~CBasedTraversal() override = default;

  /**
   * The main traversal of the CTraversal.
   * @tparam LoopBody type of the loop body
   * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
   * (x,y,z). If you need additional input from outside, please use captures (by reference).
   * @param end 3D index until interactions are processed (exclusive)
   * @param stride dimension of stride (depends on coloring)
   * @param offset initial offset
   */
  template <typename LoopBody>
  inline void cTraversal(LoopBody&& loopBody, const std::array<unsigned long, 3>& end,
                         const std::array<unsigned long, 3>& stride,
                         const std::array<unsigned long, 3>& offset = {0ul, 0ul, 0ul});

  /**
   * cutoff radius.
   */
  const ParticleFloatType _cutoff;

  /**
   * cell length in CellBlock3D.
   */
  const std::array<ParticleFloatType, 3> _cellLength;

  /**
   * overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  std::array<unsigned long, 3> _overlap;
};

template <class ParticleCell>
template <typename LoopBody>
inline void CBasedTraversal<ParticleCell>::cTraversal(LoopBody&& loopBody, const std::array<unsigned long, 3>& end,
                                                      const std::array<unsigned long, 3>& stride,
                                                      const std::array<unsigned long, 3>& offset) {
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
    const unsigned long numColors = stride[0] * stride[1] * stride[2];
    for (unsigned long col = 0; col < numColors; ++col) {
      std::array<unsigned long, 3> startWithoutOffset(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));
      std::array<unsigned long, 3> start(ArrayMath::add(startWithoutOffset, offset));

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
