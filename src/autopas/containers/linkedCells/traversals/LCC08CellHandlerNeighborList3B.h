/**
 * @file LCC08CellHandlerNeighborList3B.h
 * @author Alexander Haberl
 * @date 2024/04/15
 */

#pragma once

#include "autopas/baseFunctors/CellFunctorNeighborListBuild3B.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step for 3-Body neighbor list building.
 *
 * The base step processBaseCell() computes one set of pairwise interactions
 * between two cells for each spatial direction based on the baseIndex.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class LCC08CellHandlerNeighborList3B {
 public:
  /**
   * Constructor of the LCC08CellHandlerNeighborList3B.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC08CellHandlerNeighborList3B(PairwiseFunctor *pairwiseFunctor,
                                          const std::array<unsigned long, 3> &cellsPerDimension,
                                          const double interactionLength, const std::array<double, 3> &cellLength,
                                          const std::array<unsigned long, 3> &overlap, DataLayoutOption dataLayout,
                                          bool useNewton3)
      : _cellFunctor(pairwiseFunctor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/,
                     dataLayout, useNewton3),
        _cellPairOffsets{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlap(overlap),
        _dataLayout(dataLayout),
        _useNewton3(useNewton3) {
    computeOffsets(cellsPerDimension);
  }

  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner (=base index) of a 2x2x2 block of cells.
   * @param cells vector of all cells.
   * @param baseIndex Index respective to which box is constructed.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long baseIndex);

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) { _cellFunctor.setSortingThreshold(sortingThreshold); }

 protected:
  /**
   * Pair sets for processBaseCell().
   * Values are: offset of first cell, offset of second cell, sorting direction.
   */
  std::vector<std::tuple<unsigned long, unsigned long, std::array<double, 3>>> _cellPairOffsets;

  /**
   * Computes pairs for the block used in processBaseCell().
   * The algorithm used to generate the cell pairs can be visualized with a python script, which can be found in
   * docs/C08TraversalScheme.py
   * @param cellsPerDimension
   */
  void computeOffsets(const std::array<unsigned long, 3> &cellsPerDimension);

  /**
   * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  const std::array<unsigned long, 3> _overlap;

  /**
   * The datalayout to be used.
   */
  DataLayoutOption _dataLayout;

  /**
   * If newton3 should be used or not.
   */
  bool _useNewton3;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctorNeighborListBuild3B<ParticleCell, PairwiseFunctor, /*bidirectional*/ true> _cellFunctor;

  /**
   * Interaction length (cutoff + skin).
   */
  const double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;
};

template <class ParticleCell, class PairwiseFunctor>
inline void LCC08CellHandlerNeighborList3B<ParticleCell, PairwiseFunctor>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long baseIndex) {
  for (auto const &[offset1, offset2, r] : _cellPairOffsets) {
    const unsigned long cellIndex1 = baseIndex + offset1;
    const unsigned long cellIndex2 = baseIndex + offset2;

    if (cellIndex1 >= cells.size() or cellIndex2 >= cells.size()) {
      // check that index is not outOfBounds because we call processBaseCell on outer-most Halo-Cells as well
      continue;
    }

    ParticleCell &cell1 = cells[cellIndex1];
    ParticleCell &cell2 = cells[cellIndex2];

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
    } else {
      this->_cellFunctor.processCellPair(cell1, cell2, r);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor>
inline void LCC08CellHandlerNeighborList3B<ParticleCell, PairwiseFunctor>::computeOffsets(
    const std::array<unsigned long, 3> &cellsPerDimension) {
  using namespace autopas::utils::ArrayMath::literals;
  using std::make_pair;

  //////////////////////////////
  // @TODO: Replace following lines with vector to support asymmetric cells
  const unsigned long ov1 = _overlap[0] + 1;
  const unsigned long ov1_squared = ov1 * ov1;
  //////////////////////////////

  const std::array<unsigned long, 3> overlap_1 = _overlap + 1ul;

  std::vector<unsigned long> cellOffsets;
  cellOffsets.reserve(overlap_1[0] * overlap_1[1] * overlap_1[2]);

  _cellPairOffsets.clear();

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  // constants 0.0 and 1.0 with correct type
  const double zero = 0.0;
  const double one = 1.0;

  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
      for (unsigned long z = 0ul; z <= _overlap[2]; ++z) {
        cellOffsets.push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension));
      }
    }
  }
  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
      for (unsigned long z = 0ul; z <= _overlap[2]; ++z) {
        const unsigned long offset = cellOffsets[ov1_squared * x + ov1 * y];
        const std::array<double, 3> offsetVec =
            utils::ArrayUtils::static_cast_copy_array<double>(
                utils::ThreeDimensionalMapping::oneToThreeD(offset, cellsPerDimension)) *
            this->_cellLength;
        // origin
        {
          // check whether cell is within interaction length. distVec is the direction between borders of cells.
          const auto distVec =
              std::array<double, 3>{std::max(zero, x - one), std::max(zero, y - one), std::max(zero, z - one)} *
              _cellLength;
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            const auto baseCellVec = utils::ArrayUtils::static_cast_copy_array<double>(
                utils::ThreeDimensionalMapping::oneToThreeD(cellOffsets[z], cellsPerDimension));

            // Calculate the sorting direction from the base cell to the other cell.
            auto sortingDir = offsetVec - baseCellVec;
            if (x == 0 and y == 0 and z == 0) {
              sortingDir = {1., 1., 1.};
            }

            _cellPairOffsets.push_back(
                std::make_tuple(cellOffsets[z], offset, utils::ArrayMath::normalize(sortingDir)));
          }
        }
        // back left
        if (y != _overlap[1] and z != 0) {
          // check whether cell is within interaction length
          const auto distVec = std::array<double, 3>{std::max(zero, x - one), std::max(zero, _overlap[1] - y - one),
                                                     std::max(zero, z - one)} *
                               _cellLength;
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            const auto baseCellVec = utils::ArrayUtils::static_cast_copy_array<double>(
                utils::ThreeDimensionalMapping::oneToThreeD(cellOffsets[ov1_squared - ov1 + z], cellsPerDimension));

            // Calculate the sorting direction from the base cell to the other cell.
            auto sortingDir = offsetVec - baseCellVec;
            if (sortingDir[0] == 0 and sortingDir[1] == 0 and sortingDir[2] == 0) {
              sortingDir = {1., 1., 1.};
            }
            _cellPairOffsets.emplace_back(cellOffsets[ov1_squared - ov1 + z], offset,
                                          utils::ArrayMath::normalize(sortingDir));
          }
        }
        // front right
        if (x != _overlap[0] and (y != 0 or z != 0)) {
          // check whether cell is within interaction length
          const auto distVec = std::array<double, 3>{std::max(zero, _overlap[0] - x - one), std::max(zero, y - one),
                                                     std::max(zero, z - one)} *
                               _cellLength;
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            const auto baseCellVec =
                utils::ArrayUtils::static_cast_copy_array<double>(utils::ThreeDimensionalMapping::oneToThreeD(
                    cellOffsets[ov1_squared * _overlap[0] + z], cellsPerDimension));

            // Calculate the sorting direction from the base cell to the other cell.
            auto sortingDir = offsetVec - baseCellVec;
            if (sortingDir[0] == 0 and sortingDir[1] == 0 and sortingDir[2] == 0) {
              sortingDir = {1., 1., 1.};
            }
            sortingDir = utils::ArrayMath::normalize(sortingDir);
            _cellPairOffsets.emplace_back(cellOffsets[ov1_squared * _overlap[0] + z], offset,
                                          utils::ArrayMath::normalize(sortingDir));
          }
        }
        // back right
        if (y != _overlap[1] and x != _overlap[0] and z != 0) {
          // check whether cell is within interaction length
          const auto distVec = std::array<double, 3>{std::max(zero, _overlap[0] - x - one),
                                                     std::max(zero, _overlap[1] - y - one), std::max(zero, z - one)} *
                               _cellLength;
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            const auto baseCellVec =
                utils::ArrayUtils::static_cast_copy_array<double>(utils::ThreeDimensionalMapping::oneToThreeD(
                    cellOffsets[ov1_squared * ov1 - ov1 + z], cellsPerDimension));

            // Calculate the sorting direction from the base cell to the other cell.
            auto sortingDir = offsetVec - baseCellVec;
            if (sortingDir[0] == 0 and sortingDir[1] == 0 and sortingDir[2] == 0) {
              sortingDir = {1., 1., 1.};
            }

            _cellPairOffsets.emplace_back(cellOffsets[ov1_squared * ov1 - ov1 + z], offset,
                                          utils::ArrayMath::normalize(sortingDir));
          }
        }
      }
    }
  }
}

}  // namespace autopas
