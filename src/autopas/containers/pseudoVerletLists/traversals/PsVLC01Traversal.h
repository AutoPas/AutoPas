/**
 * @file PsVLC01Traversal.h
 * @author Lars Doll
 * @date 12.01.2026
 */

#pragma once

#include "PsVLTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctorPsVL.h"
#include "autopas/containers/cellTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"

namespace autopas {

/**
 * This class provides the psVl_c01 traversal.
 * @tparam ParticleCell the type of cells.
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class PsVLC01Traversal : public C01BasedTraversal<ParticleCell, PairwiseFunctor, 3>,
                         public PsVLTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param pairwiseFunctor The functor that defines the interaction of particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length in CellBlock3D
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * in that case the interactionLength is needed!
   */
  explicit PsVLC01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                            const double interactionLength, const std::array<double, 3> &cellLength,
                            DataLayoutOption dataLayout, bool useNewton3)
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, 3>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                            dataLayout, useNewton3),
        _cellFunctor(pairwiseFunctor, interactionLength, dataLayout, useNewton3) {
    computeOffsets();
  }

  /**
   * Computes pairs used in processBaseCell().
   */
  void computeOffsets();

  void traverseParticles() override;

  /**
   * C01 traversals are only usable if useNewton3 is disabled.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override {
    if (this->_useNewton3 == false && this->_dataLayout == DataLayoutOption::aos) {
      return true;
    }
    return false;
  }
  /**
   * Getter.
   * @return
   */
  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::psvl_c01; }

  /**
   * Sets the orientationList.
   * @param list
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell>>> &list) override;

  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  using CellOffsetsType = std::vector<std::vector<std::pair<long, std::array<double, 3>>>>;

  /**
   * Computes all interactions between the base
   * cell and adjacent cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  inline void processBaseCell(unsigned long x, unsigned long y, unsigned long z);

  /**
   * Pairwise implementation of processBaseCell().
   * @copydoc processBaseCell()
   */
  inline void processBaseCellPairwise(unsigned long x, unsigned long y, unsigned long z);

  /**
   * Pairwise implementation of computeOffsets().
   * @copydoc computeOffsets()
   */
  void computePairwiseOffsets();

  /**
   * Pairs or triplets for processBaseCell().
   * @note std::map not applicable since ordering arising from insertion is important for later processing!
   */
  CellOffsetsType _cellOffsets;

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctorPsVL<ParticleCell, PairwiseFunctor,
                            /*bidirectional*/ false>
      _cellFunctor;
};

template <class ParticleCell, class PairwiseFunctor>
void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell>>> &list) {
  PsVLTraversalInterface<ParticleCell>::setOrientationList(list);
  _cellFunctor.setOrientationList(list);
}

template <class ParticleCell, class PairwiseFunctor>
inline void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::computeOffsets() {
  computePairwiseOffsets();
}

template <class ParticleCell, class PairwiseFunctor>
inline void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::computePairwiseOffsets() {
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
inline void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::processBaseCell(unsigned long x, unsigned long y,
                                                                             unsigned long z) {
  processBaseCellPairwise(x, y, z);
}

template <class ParticleCell, class PairwiseFunctor>
inline void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::processBaseCellPairwise(unsigned long x, unsigned long y,
                                                                                     unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  for (const auto &slice : _cellOffsets) {
    for (auto const &[offset, r] : slice) {
      const unsigned long otherIndex = baseIndex + offset;

      if (baseIndex == otherIndex) {
        this->_cellFunctor.processCell(baseIndex);
      } else {
        this->_cellFunctor.processCellPair(baseIndex, otherIndex, r);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor>
inline void PsVLC01Traversal<ParticleCell, PairwiseFunctor>::traverseParticles() {
  auto &cells = *(this->_cells);
  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(x, y, z); });
}
}  // namespace autopas
