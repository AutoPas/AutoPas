/**
 * @file C01BasedTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <algorithm>
#include <array>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "ColorBasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/checkFunctorType.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c01 base step.
 *
 * The traversal is defined in the function c01Traversal and uses 1 color. Interactions between cells are allowed
 * only if particles of only one cell are modified. This means that newton3 optimizations are NOT allowed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of particles.
 * @tparam dataLayout indicates usage of SoA
 */
template <class ParticleCell, class Functor, int collapseDepth = 3>
class C01BasedTraversal : public ColorBasedTraversal<ParticleCell, Functor, collapseDepth> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit C01BasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor, double interactionLength,
                             const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : ColorBasedTraversal<ParticleCell, Functor, collapseDepth>(dims, functor, interactionLength, cellLength,
                                                                  dataLayout, useNewton3) {}

  /**
   * Computes all combinations of cells used in processBaseCell()
   */
  void computeOffsets();

 protected:
  /**
   * The main traversal of the C01Traversal.
   * This provides the structure of the loops and its parallelization.
   * @tparam LoopBody
   * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
   * (x,y,z). If you need additional input from outside, please use captures (by reference).
   */
  template <typename LoopBody>
  inline void c01Traversal(LoopBody &&loopBody);

  /**
   * Pairwise implementation of computeOffsets().
   * @copydoc computeOffsets()
   */
  void computePairwiseOffsets();

  /**
   * Triwise implementation of computeOffsets().
   * @copydoc computeOffsets()
   */
  void computeTriwiseOffsets();

  // CelllOffsets type for either Pairwise or Triwise Functors.
  using CellOffsetsType = std::conditional_t<decltype(utils::isPairwiseFunctor<Functor>())::value,
                                             std::vector<std::vector<std::pair<long, std::array<double, 3>>>>,
                                             std::vector<std::tuple<long, long, std::array<double, 3>>>>;

  /**
   * Pairs or triplets for processBaseCell().
   * @note std::map not applicable since ordering arising from insertion is important for later processing!
   */
  CellOffsetsType _cellOffsets;
};

template <class ParticleCell, class Functor, int collapseDepth>
template <typename LoopBody>
inline void C01BasedTraversal<ParticleCell, Functor, collapseDepth>::c01Traversal(LoopBody &&loopBody) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto offset = this->_overlap;
  const auto end = this->_cellsPerDimension - this->_overlap;
  this->colorTraversal(std::forward<LoopBody>(loopBody), end, {1ul, 1ul, 1ul}, offset);
}

template <class ParticleCell, class Functor, int collapseDepth>
void C01BasedTraversal<ParticleCell, Functor, collapseDepth>::computeOffsets() {
  if constexpr (utils::isPairwiseFunctor<Functor>()) {
    computePairwiseOffsets();
  } else if (utils::isTriwiseFunctor<Functor>()) {
    computeTriwiseOffsets();
  } else {
    utils::ExceptionHandler::exception("C01BasedTraversal::computeOffsets(): Functor is not valid.");
  }
}

template <class ParticleCell, class Functor, int collapseDepth>
void C01BasedTraversal<ParticleCell, Functor, collapseDepth>::computePairwiseOffsets() {
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

template <class ParticleCell, class Functor, int collapseDepth>
void C01BasedTraversal<ParticleCell, Functor, collapseDepth>::computeTriwiseOffsets() {
  using namespace utils::ArrayMath::literals;
  // Reserve approximately. Overestimates more for larger overlap.
  const int cubeSize = this->_overlap[0] * this->_overlap[1] * this->_overlap[2];
  _cellOffsets.reserve(cubeSize * cubeSize / 4);

  // Helper function to get minimal distance between two cells
  auto cellDistance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
    return std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * this->_cellLength[0],
                                 std::max(0l, (std::abs(y1 - y2) - 1l)) * this->_cellLength[1],
                                 std::max(0l, (std::abs(z1 - z2) - 1l)) * this->_cellLength[2]};
  };

  const auto interactionLengthSquare{this->_interactionLength * this->_interactionLength};
  _cellOffsets.emplace_back(0, 0, std::array<double, 3>{1., 1., 1.});

  // offsets for the first cell
  for (long x1 = -this->_overlap[0]; x1 <= static_cast<long>(this->_overlap[0]); ++x1) {
    for (long y1 = -this->_overlap[1]; y1 <= static_cast<long>(this->_overlap[1]); ++y1) {
      for (long z1 = -this->_overlap[2]; z1 <= static_cast<long>(this->_overlap[2]); ++z1) {
        // check distance between base cell and cell 1
        const auto dist01 = cellDistance(0l, 0l, 0l, x1, y1, z1);

        const double distSquare = utils::ArrayMath::dot(dist01, dist01);
        if (distSquare > interactionLengthSquare) continue;

        // offsets for cell 2
        for (long x2 = -this->_overlap[0]; x2 <= static_cast<long>(this->_overlap[0]); ++x2) {
          for (long y2 = -this->_overlap[1]; y2 <= static_cast<long>(this->_overlap[1]); ++y2) {
            for (long z2 = -this->_overlap[2]; z2 <= static_cast<long>(this->_overlap[2]); ++z2) {
              // check distance between cell 1 and cell 2
              const auto dist12 = cellDistance(x1, y1, z1, x2, y2, z2);

              const double dist12Squared = utils::ArrayMath::dot(dist12, dist12);
              if (dist12Squared > interactionLengthSquare) continue;

              // check distance between base cell and cell 2
              const auto dist02 = cellDistance(0l, 0l, 0l, x2, y2, z2);

              const double dist02Squared = utils::ArrayMath::dot(dist02, dist02);
              if (dist02Squared > interactionLengthSquare) continue;

              const long offset1 = utils::ThreeDimensionalMapping::threeToOneD(
                  x1, y1, z1, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));

              const long offset2 = utils::ThreeDimensionalMapping::threeToOneD(
                  x2, y2, z2, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));

              // Only add unique combinations. E.g.: (5, 8) == (8, 5)
              if (offset2 <= offset1) continue;

              // sorting direction from base cell to the first different cell
              std::array<double, 3> sortDirection{};
              if (offset1 == 0) {
                sortDirection = {x2 * this->_cellLength[0], y2 * this->_cellLength[1], z2 * this->_cellLength[2]};
              } else {
                sortDirection = {x1 * this->_cellLength[0], y1 * this->_cellLength[1], z1 * this->_cellLength[2]};
              }
              _cellOffsets.emplace_back(offset1, offset2, utils::ArrayMath::normalize(sortDirection));
            }
          }
        }
      }
    }
  }
}

}  // namespace autopas
