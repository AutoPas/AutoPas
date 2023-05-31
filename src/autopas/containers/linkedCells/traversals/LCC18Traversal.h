/**
 * @file LCC18Traversal.h
 * @author nguyen
 * @date 06.09.2018
 */

#pragma once

#include "LCTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C18BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the lc_c18 traversal.
 *
 * The traversal uses the c18 base step performed on every single cell.
 * \image html C18.png "C18 base step in 2D. (dark blue cell = base cell)"
 * Since these steps overlap a domain coloring with eighteen colors is applied.
 * \image html C18_domain.png "C18 domain coloring in 2D. 6 colors are required."
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC18Traversal : public C18BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                       public LCTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the lc_c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC18Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                          const double interactionLength, const std::array<double, 3> &cellLength)
      : C18BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _cellFunctor(pairwiseFunctor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/) {
    computeOffsets();
  }

  void traverseParticlePairs() override;

  /**
   * Computes all interactions between the base
   * cell and adjacent cells with greater a ID.
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c18; }

  /**
   * C18 traversal is always usable.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override { return true; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

 private:
  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, dataLayout, useNewton3>
      _cellFunctor;

  /**
   * Type of an array containing offsets relative to the base cell and correspondent normalized 3d relationship vectors.
   * The vectors (aka std::array<double,3>) describe the imaginative line connecting the center of the base cell and the
   * center of the cell defined by the offset. It is used for sorting.
   */
  using offsetArray_t = std::vector<std::pair<unsigned long, std::array<double, 3>>>;

  /**
   * Pairs for processBaseCell(). overlap[0] x overlap[1] offsetArray_t for each special case in x and y direction.
   */
  std::vector<std::vector<offsetArray_t>> _cellOffsets;

  /**
   * Returns the index in the offsets array for the given position.
   * @param pos current position in dimension dim
   * @param dim current dimension
   * @return Index for the _cellOffsets Array.
   */
  unsigned long getIndex(const unsigned long pos, const unsigned int dim) const;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::computeOffsets() {
  _cellOffsets.resize(2 * this->_overlap[1] + 1, std::vector<offsetArray_t>(2 * this->_overlap[0] + 1));
  const std::array<long, 3> _overlap_s = utils::ArrayUtils::static_cast_array<long>(this->_overlap);

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  for (long z = 0l; z <= _overlap_s[2]; ++z) {
    for (long y = -_overlap_s[1]; y <= _overlap_s[1]; ++y) {
      for (long x = -_overlap_s[0]; x <= _overlap_s[0]; ++x) {
        const long offset = utils::ThreeDimensionalMapping::threeToOneD(
            x, y, z, utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension));

        if (offset < 0l) {
          continue;
        }
        // add to each applicable special case
        for (long yArray = -_overlap_s[1]; yArray <= _overlap_s[1]; ++yArray) {
          if (std::abs(yArray + y) <= _overlap_s[1]) {
            for (long xArray = -_overlap_s[0]; xArray <= _overlap_s[0]; ++xArray) {
              if (std::abs(xArray + x) <= _overlap_s[0]) {
                std::array<double, 3> pos = {};
                pos[0] = std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0];
                pos[1] = std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1];
                pos[2] = std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2];
                // calculate distance between base cell and other cell
                const double distSquare = utils::ArrayMath::dot(pos, pos);
                // only add cell offset if cell is within cutoff radius
                if (distSquare <= interactionLengthSquare) {
                  _cellOffsets[yArray + _overlap_s[1]][xArray + _overlap_s[0]].push_back(
                      std::make_pair(offset, utils::ArrayMath::normalize(pos)));
                }
              }
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
unsigned long LCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::getIndex(
    const unsigned long pos, const unsigned int dim) const {
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

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void LCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  const unsigned long xArray = getIndex(x, 0);
  const unsigned long yArray = getIndex(y, 1);

  ParticleCell &baseCell = cells[baseIndex];

  offsetArray_t &offsets = this->_cellOffsets[yArray][xArray];
  for (auto const &[offset, r] : offsets) {
    unsigned long otherIndex = baseIndex + offset;
    ParticleCell &otherCell = cells[otherIndex];

    if (baseIndex == otherIndex) {
      if (!this->_onlyDirty || baseCell.getDirty()) {
        this->_cellFunctor.processCell(baseCell);
      }
    } else {
      if (!this->_onlyDirty || baseCell.getDirty() || otherCell.getDirty()) {
        this->_cellFunctor.processCellPair(baseCell, otherCell, r);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC18Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  this->c18Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
