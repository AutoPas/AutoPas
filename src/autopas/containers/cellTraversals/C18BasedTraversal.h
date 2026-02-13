/**
 * @file C18BasedTraversal.h
 * @author nguyen
 * @date 06.09.18
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

#include "autopas/containers/cellTraversals/ColorBasedTraversal.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c18 base step.
 * The traversal is defined in the function c18Traversal and uses 18 colors, such that interactions between the base
 * cell and all adjacent cells with greater ID are safe, even when using newton3 optimizations.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class Functor>
class C18BasedTraversal : public ColorBasedTraversal<ParticleCell, Functor> {
 public:
  /**
   * Constructor of the lc_c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit C18BasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor, const double interactionLength,
                             const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : ColorBasedTraversal<ParticleCell, Functor>(dims, functor, interactionLength, cellLength, dataLayout,
                                                   useNewton3) {}

 protected:
  /**
   * The main traversal of the C18Traversal.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   *
   * @tparam allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the c18 step if allCells is false, iteration will not occur over the last layer of
   * cells (for overlap=1) (in x, y and z direction).
   */
  template <bool allCells, typename LoopBody>
  inline void c18Traversal(LoopBody &&loopBody);

  void computeOffsets();

  using offsetArray_t = std::vector<std::pair<unsigned long, std::array<double, 3>>>;
  std::vector<std::vector<offsetArray_t>> _cellOffsets;

  [[nodiscard]] unsigned long getIndex(unsigned long pos, unsigned int dim) const;
};

template <class ParticleCell, class Functor>
template <bool allCells, typename LoopBody>
inline void C18BasedTraversal<ParticleCell, Functor>::c18Traversal(LoopBody &&loopBody) {
  const std::array<unsigned long, 3> stride = {
      2ul * this->_overlap[0] + 1ul,
      2ul * this->_overlap[1] + 1ul,
      this->_overlap[2] + 1ul,
  };
  auto end(this->_cellsPerDimension);
  if (not allCells) {
    end[2] -= this->_overlap[2];
  }
  this->colorTraversal(std::forward<LoopBody>(loopBody), end, stride);
}

template <class ParticleCell, class Functor>
inline void C18BasedTraversal<ParticleCell, Functor>::computeOffsets() {
  _cellOffsets.resize(2 * this->_overlap[1] + 1, std::vector<offsetArray_t>(2 * this->_overlap[0] + 1));
  const std::array<long, 3> _overlap_s = utils::ArrayUtils::static_cast_copy_array<long>(this->_overlap);

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  for (long z = 0l; z <= _overlap_s[2]; ++z) {
    for (long y = -_overlap_s[1]; y <= _overlap_s[1]; ++y) {
      for (long x = -_overlap_s[0]; x <= _overlap_s[0]; ++x) {
        const long offset = utils::ThreeDimensionalMapping::threeToOneD(
            x, y, z, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));

        if (offset < 0l) {
          continue;
        }
        // add to each applicable special case
        for (long yArray = -_overlap_s[1]; yArray <= _overlap_s[1]; ++yArray) {
          if (std::abs(yArray + y) <= _overlap_s[1]) {
            for (long xArray = -_overlap_s[0]; xArray <= _overlap_s[0]; ++xArray) {
              if (std::abs(xArray + x) <= _overlap_s[0]) {
                const std::array<double, 3> pos = {
                    std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0],
                    std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1],
                    std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2],
                };
                // calculate distance between the borders of the base cell and the other cell
                const double distSquare = utils::ArrayMath::dot(pos, pos);
                // only add cell offset if cell is within cutoff radius
                if (distSquare <= interactionLengthSquare) {
                  // Calculate the sorting direction from the base cell and the other cell by use of the offset (x, y,
                  // z).
                  // Note: We have to calculate the sorting direction separately from pos, since pos is the offset
                  // between the borders of cells. For neighbouring cells this would be 0 and the sorting direction
                  // would be wrong.
                  std::array<double, 3> sortingDir = {static_cast<double>(x) * this->_cellLength[0],
                                                      static_cast<double>(y) * this->_cellLength[1],
                                                      static_cast<double>(z) * this->_cellLength[2]};
                  if (x == 0 and y == 0 and z == 0) {
                    sortingDir = {1., 1., 1.};
                  }

                  _cellOffsets[yArray + _overlap_s[1]][xArray + _overlap_s[0]].push_back(
                      std::make_pair(offset, utils::ArrayMath::normalize(sortingDir)));
                }
              }
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class Functor>
unsigned long C18BasedTraversal<ParticleCell, Functor>::getIndex(const unsigned long pos,
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

}  // namespace autopas
