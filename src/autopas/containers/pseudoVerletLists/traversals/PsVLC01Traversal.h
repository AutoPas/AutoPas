/**
 * @file PsVLC01Traversal.h
 * @author Lars Doll
 * @date 01.01.2026
 */

#pragma once

#include "PsVLTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/baseFunctors/CellFunctorPsVL.h"
#include "autopas/containers/cellTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/checkFunctorType.h"

namespace autopas {

/**
 * This class provides the psVl_c01 traversal.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class PsVLC01Traversal : public C01BasedTraversal<ParticleCell, PairwiseFunctor, 3>,
                         public PsVLTraversalInterface<ParticleCell> {
 public:
  explicit PsVLC01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                            const double interactionLength, const std::array<double, 3> &cellLength,
                            DataLayoutOption dataLayout, bool useNewton3)
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, 3>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                            dataLayout, useNewton3),
        _cellFunctor(pairwiseFunctor, interactionLength, dataLayout, useNewton3) {
    computeOffsets();
  }

 private:
  /**
   * Computes pairs used in processBaseCell().
   */
  void computeOffsets();

  /**
   * Returns the index in the offset array for the given position.
   * @param pos current position in dimension dim.
   * @param dim current dimension.
   * @return Index for the _cellOffsets Array.
   */
  [[nodiscard]] unsigned long getIndex(unsigned long pos, unsigned int dim) const;

  /**
   * Getter.
   */
  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::psvl_c01; }

  /**
   * C01 traversal is always usable.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override { return true; }

  /**
   * Sets the orientationList.
   * @param list
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell>>> &list) override;

  void setSortingThreshold(size_t sortingThreshold) override {}

  /**
   * Computes all interactions between the base cell and adjacent cells.
   * @param x X-index of the base cell.
   * @param y Y-index of the base cell.
   * @param z Z-index of the base cell.
   */
  void processBaseCell(unsigned long x, unsigned long y, unsigned long z);

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctorPsVL<ParticleCell, PairwiseFunctor,
                            /*bidirectional*/ true>
      _cellFunctor;

  /**
   * Type of array containing offsets relative to the base cell and correspondent normalized 3d relationship vectors.
   * The vectors (aka std::array<double,3>) describe the imaginative line connecting the center of the base cell and the
   * center of the cell defined by the offset. It is used for sorting.
   */
  using offsetArray_t = std::vector<std::pair<unsigned long, std::array<double, 3>>>;

  /**
   * Pairs for processBaseCell(). overlap[0] x overlap[1] offsetArray_t for each special case in x and y direction.
   */
  std::vector<offsetArray_t> _cellOffsets;
};

template <class ParticleCell, class PairwiseFunctor>
void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell>>> &list) {
  PsVLTraversalInterface<ParticleCell>::setOrientationList(list);
  _cellFunctor.setOrientationList(list);
}

template <class ParticleCell, class Functor>
inline void PsVLC01Traversal<ParticleCell, Functor>::computeOffsets() {
  _cellOffsets.resize(2 * this->_overlap[0] + 1);

  const auto interactionLengthSquare{this->_interactionLength * this->_interactionLength};

  for (long x = -this->_overlap[0]; x <= 0l; ++x) {
    for (long y = -this->_overlap[1]; y <= static_cast<long>(this->_overlap[1]); ++y) {
      for (long z = -this->_overlap[2]; z <= static_cast<long>(this->_overlap[2]); ++z) {
        const std::array<double, 3> pos = {
            std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0],
            std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1],
            std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2],
        };
        const double distSquare = utils::ArrayMath::dot(pos, pos);
        if (distSquare <= interactionLengthSquare) {
          const long currentOffset = utils::ThreeDimensionalMapping::threeToOneD(
              x, y, z, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));
          const bool containCurrentOffset =
              std::any_of(_cellOffsets[x + this->_overlap[0]].cbegin(), _cellOffsets[x + this->_overlap[0]].cend(),
                          [currentOffset](const auto &e) { return e.first == currentOffset; });
          if (containCurrentOffset) {
            continue;
          }
          for (long ix = x; ix <= std::abs(x); ++ix) {
            const long offset = utils::ThreeDimensionalMapping::threeToOneD(
                ix, y, z, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));
            const size_t index = ix + this->_overlap[0];

            // Calculate the sorting direction from the base cell (x, y, z) and the other cell by use of the offset (ix,
            // y, z).
            std::array<double, 3> sortingDir = {static_cast<double>(ix) * this->_cellLength[0],
                                                static_cast<double>(y) * this->_cellLength[1],
                                                static_cast<double>(z) * this->_cellLength[2]};

            // the offset to the current cell itself is zero.
            if (ix == 0 and y == 0 and z == 0) {
              sortingDir = {1., 1., 1.};
            }
            sortingDir = utils::ArrayMath::normalize(sortingDir);

            if (y == 0l and z == 0l) {
              // make sure center of slice is always at the beginning
              _cellOffsets[index].insert(_cellOffsets[index].cbegin(), std::make_pair(offset, sortingDir));
            } else {
              _cellOffsets[index].emplace_back(offset, sortingDir);
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor>
unsigned long PsVLC01Traversal<ParticleCell, PairwiseFunctor>::getIndex(const unsigned long pos,
                                                                        const unsigned int dim) const {
  unsigned long index;
  if (pos < this->_overlap[dim]) {
    index = pos;
  } else if (pos < this->_cellsPerDimension[dim] - this->_overlap[dim]) {
    index = this->_overlap[dim];
  } else {
    index = pos - this->_cellsPerDimension[dim] + 2 * this->_overlap[dim] + 1ul;
  }
  return index;
}

template <class ParticleCell, class PairwiseFunctor>
void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::processBaseCell(unsigned long x, unsigned long y,
                                                                      unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  const unsigned long xArray = getIndex(x, 0);

  offsetArray_t &offsets = _cellOffsets[xArray];

  for (auto const &[offset, r] : offsets) {
    const unsigned long otherIndex = baseIndex + offset;

    if (baseIndex == otherIndex) {
      this->_cellFunctor.processCell(baseIndex);
    } else {
      this->_cellFunctor.processCellPair(baseIndex, otherIndex, r);
    }
  }
}

}  // namespace autopas