/**
 * @file CBasedTraversal.h
 * @author C. Menges
 * @date 26.04.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the base for traversals using base steps based on cell coloring.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam collapseDepth Set the depth of loop collapsion for OpenMP. Loop variables from outer to inner loop: z,y,x
 */
template <class ParticleCell, class PairwiseFunctor, int collapseDepth = 3>
class CBasedTraversal : public CellPairTraversal<ParticleCell> {
 protected:
  /**
   * Constructor of the CBasedTraversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit CBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                           double interactionLength, const std::array<double, 3> &cellLength,
                           DataLayoutOption::Value dataLayout, bool useNewton3)
      : CellPairTraversal<ParticleCell>(dims, dataLayout, useNewton3),
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _dataLayoutConverter(pairwiseFunctor, dataLayout) {
    for (unsigned int d = 0; d < 3; d++) {
      _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
    }
  }

  /**
   * Destructor of CBasedTraversal.
   */
  ~CBasedTraversal() override = default;

 public:
  /**
   * load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  void initTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
      /// @todo find a condition on when to use omp or when it is just overhead
      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.loadDataLayout(cells[i]);
      }
    }
  }

  /**
   * write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
      /// @todo find a condition on when to use omp or when it is just overhead
      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.storeDataLayout(cells[i]);
      }
    }
  }

 protected:
  /**
   * The main traversal of the CTraversal.
   * @tparam LoopBody type of the loop body
   * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
   * (x,y,z). If you need additional input from outside, please use captures (by reference).
   * @param end 3D index until interactions are processed (exclusive).
   * @param stride Distance (in cells) to the next cell of the same color.
   * @param offset initial offset (in cells) in which cell to start the traversal.
   */
  template <typename LoopBody>
  inline void cTraversal(LoopBody &&loopBody, const std::array<unsigned long, 3> &end,
                         const std::array<unsigned long, 3> &stride,
                         const std::array<unsigned long, 3> &offset = {0ul, 0ul, 0ul});

  /**
   * This method is called when the color during the traversal has changed.
   *
   * @param newColor The new current color.
   */
  virtual void notifyColorChange(unsigned long newColor){};

  /**
   * Interaction length (cutoff + skin).
   */
  const double _interactionLength;

  /**
   * cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;

  /**
   * overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  std::array<unsigned long, 3> _overlap;

 private:
  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, int collapseDepth>
template <typename LoopBody>
inline void CBasedTraversal<ParticleCell, PairwiseFunctor, collapseDepth>::cTraversal(
    LoopBody &&loopBody, const std::array<unsigned long, 3> &end, const std::array<unsigned long, 3> &stride,
    const std::array<unsigned long, 3> &offset) {
  using namespace autopas::utils::ArrayMath::literals;
  AUTOPAS_OPENMP(parallel) {
    const unsigned long numColors = stride[0] * stride[1] * stride[2];
    for (unsigned long col = 0; col < numColors; ++col) {
      AUTOPAS_OPENMP(single) {
        // barrier at omp for of previous loop iteration, so fine to change it for everyone!
        notifyColorChange(col);
        // implicit barrier at end of function.
      }
      const std::array<unsigned long, 3> startWithoutOffset(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));
      const std::array<unsigned long, 3> start(startWithoutOffset + offset);

      // intel compiler demands following:
      const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
      const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
      const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];
      if (collapseDepth == 2) {
        AUTOPAS_OPENMP(for schedule(dynamic, 1) collapse(2))
        for (unsigned long z = start_z; z < end_z; z += stride_z) {
          for (unsigned long y = start_y; y < end_y; y += stride_y) {
            for (unsigned long x = start_x; x < end_x; x += stride_x) {
              // Don't exchange order of execution (x must be last!), it would break other code
              loopBody(x, y, z);
            }
          }
        }
      } else {
        AUTOPAS_OPENMP(for schedule(dynamic, 1) collapse(3))
        for (unsigned long z = start_z; z < end_z; z += stride_z) {
          for (unsigned long y = start_y; y < end_y; y += stride_y) {
            for (unsigned long x = start_x; x < end_x; x += stride_x) {
              // Don't exchange order of execution (x must be last!), it would break other code
              loopBody(x, y, z);
            }
          }
        }
      }
    }
  }
}

}  // namespace autopas
