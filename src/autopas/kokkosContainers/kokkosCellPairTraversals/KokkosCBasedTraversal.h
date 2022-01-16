/**
 * @file KokkosCBasedTraversal.h
 * @author lgaertner
 * @date 29.12.2021
 */

#pragma once

#include "autopas/kokkosContainers/kokkosCellPairTraversals/KokkosCellPairTraversal.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using base steps based on cell coloring.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class KokkosCBasedTraversal : public KokkosCellPairTraversal<ParticleCell> {
 protected:
  /**
   * Constructor of the KokkosCBasedTraversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit KokkosCBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                 const double interactionLength, const std::array<double, 3> &cellLength)
      : KokkosCellPairTraversal<ParticleCell>(dims),
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _dataLayoutConverter(pairwiseFunctor),
        _overlap() {
    for (unsigned int d = 0; d < 3; d++) {
      _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
    }
  }

  /**
   * Destructor of KokkosCBasedTraversal.
   */
  ~KokkosCBasedTraversal() override = default;

 public:
  /**
   * load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  void initTraversal() override {
    // TODO lgaertner when SoA is supported
    //    if (this->_cells) {
    //      auto &cells = *(this->_cells);
    //      for (size_t i = 0; i < cells.size(); ++i) {
    //        _dataLayoutConverter.loadDataLayout(cells[i]);
    //      }
    //    }
  }

  /**
   * write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    // TODO lgaertner when SoA is supported
    //    if (this->_cells) {
    //      auto &cells = *(this->_cells);
    //      for (size_t i = 0; i < cells.size(); ++i) {
    //        _dataLayoutConverter.storeDataLayout(cells[i]);
    //      }
    //    }
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
  std::array<unsigned long, 3> _overlap{};

 private:
  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <typename LoopBody>
inline void KokkosCBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::cTraversal(
    LoopBody &&loopBody, const std::array<unsigned long, 3> &end, const std::array<unsigned long, 3> &stride,
    const std::array<unsigned long, 3> &offset) {
  const unsigned long numColors = stride[0] * stride[1] * stride[2];
  for (unsigned long col = 0; col < numColors; ++col) {
    //#if defined(AUTOPAS_OPENMP)
    //#pragma omp single
    //#endif
    {
      // barrier at omp for of previous loop iteration, so fine to change it for everyone!
      notifyColorChange(col);
      // implicit barrier at end of function.
    }
    std::array<unsigned long, 3> startWithoutOffset(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));
    std::array<unsigned long, 3> start(utils::ArrayMath::add(startWithoutOffset, offset));

    // intel compiler demands following:
    const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
    const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];

    std::vector<std::array<size_t, 3>> ranges = {};
    for (unsigned long z = start[2]; z <= end[2]; z += stride[2]) {
      for (unsigned long y = start[1]; y <= end[1]; y += stride[1]) {
        for (unsigned long x = start[0]; x <= end[0]; x += stride[0]) {
          ranges.push_back({x, y, z});
        }
      }
    }

    Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>> rangePolicy(0, ranges.size());

    Kokkos::parallel_for(
        rangePolicy, KOKKOS_LAMBDA(const size_t i) { loopBody(ranges[i][0], ranges[i][1], ranges[i][2]); });

    Kokkos::fence();
    //    const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];
    //      for (unsigned long z = start_z; z < end_z; z += stride_z) {
    //        for (unsigned long y = start_y; y < end_y; y += stride_y) {
    //          for (unsigned long x = start_x; x < end_x; x += stride_x) {
    //            // Don't exchange order of execution (x must be last!), it would break other code
    //            loopBody(x, y, z);
    //          }
    //        }
    //      }
  }
}

}  // namespace autopas
