/**
 * @file HGC08CellHandler.h
 * @author Alexander Glas
 * @date 12.04.2026
 */

#pragma once

#include <cmath>
#include <limits>

#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * This class provides the base for traversals using the hgc08 base step.
 *
 * The base step processBaseCell() computes one set of pairwise interactions
 * between two cells for each spatial direction based on the baseIndex. It also handles top-down the inter-level
 * interactions in the same spatial region. It only works for hierarchical grids with fitted grids.
 *
 * A version working for non-fitted grids was implemented, but deemed not worth it. It can be found at
 * https://github.com/AutoPas/AutoPas/pull/1136
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam PairwiseFunctor_T The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor_T>
class HGC08CellHandler : public LCC08CellHandler<ParticleCell_T, PairwiseFunctor_T> {
 public:
  /**
   * Constructor of the HGC08CellHandler.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @param cellBlocks The cell blocks of the hierarchical grid, ordered from lowest to highest level.
   * @param interactionLengthsSquared The squared interaction lengths from each lower level to upperLevel, ordered
   * from lowest to highest level.
   * @param upperLevel The upper level of the hierarchical grid
   * @todo Check if passing the cutoff instead of the interaction length to the cell functor works.
   */
  explicit HGC08CellHandler(PairwiseFunctor_T *pairwiseFunctor, const std::array<unsigned long, 3> &cellsPerDimension,
                            double interactionLength, const std::array<double, 3> &cellLength,
                            const std::array<unsigned long, 3> &overlap, DataLayoutOption dataLayout, bool useNewton3,
                            const std::vector<internal::CellBlock3D<ParticleCell_T> *> &cellBlocks,
                            const std::vector<double> &interactionLengthsSquared, const size_t upperLevel)
      : LCC08CellHandler<ParticleCell_T, PairwiseFunctor_T>(pairwiseFunctor, cellsPerDimension, interactionLength,
                                                            cellLength, overlap, dataLayout, useNewton3),
        _cellBlocks(cellBlocks),
        _interactionLengthsSquared(interactionLengthsSquared),
        _upperLevel(upperLevel) {
    using namespace autopas::utils::ArrayMath::literals;
    _shiftLength = _cellBlocks[0]->getCellLength() * 0.01;
    if (this->_upperLevel == 0) {
      return;
    }
    _cellPairOffsetsPerLevel.resize(_upperLevel);
    for (auto const &[offset1, offset2, r] : this->_cellPairOffsets) {
      computePairwiseCellOffsetsHGC08(offset1, offset2);
      if (offset1 != offset2) {
        computePairwiseCellOffsetsHGC08(offset2, offset1);
      }
    }
  }

  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner (=base index) of a 2x2x2 block of cells. Also computes top-down inter-level interactions in that
   * spatial region.
   * @param baseIndex Index of the upper-level cell respective to which box is constructed.
   */
  void processBaseCell(size_t baseIndex);

 protected:
  /**
   * Particle type handled by the particle cells.
   */
  using Particle = typename ParticleCell_T::ParticleType;
  /**
   * Cell block type used by the hierarchical-grid levels.
   */
  using CellBlock = internal::CellBlock3D<ParticleCell_T>;
  /**
   * The upper level of the hierarchical grid.
   */
  const size_t _upperLevel;
  /**
   * The cell blocks of the hierarchical grid, ordered from lowest to highest level.
   */
  std::vector<internal::CellBlock3D<ParticleCell_T> *> _cellBlocks;
  /**
   * The squared interaction lengths from each lower level to upperLevel, ordered from lowest to highest level.
   */
  std::vector<double> _interactionLengthsSquared;
  /**
   * Shift length used to shift the boundaries for determining which lower-level cells are within each upper-level cell
   * to prevent spilling into neighboring cells due to floating point errors. Must be smaller than the smallest cell
   * length to prevent skipping cells.
   */
  std::array<double, 3> _shiftLength;

  std::vector<std::vector<LCC08CellHandlerUtility::OffsetPairSorting>> _cellPairOffsetsPerLevel{};

  void computePairwiseCellOffsetsHGC08(const unsigned long offset1, const unsigned long offset2) {
    using namespace autopas::utils::ArrayMath::literals;
    // Get bounds of  cell2. Use shifted corners so exact boundary points do not spill into the neighboring cell.
    auto [lowCornerCell2, highCornerCell2] = _cellBlocks[_upperLevel]->getCellBoundingBox(offset2);
    highCornerCell2 = {highCornerCell2[0] - _shiftLength[0], highCornerCell2[1] - _shiftLength[1],
                       highCornerCell2[2] - _shiftLength[2]};

    lowCornerCell2 = {lowCornerCell2[0] + _shiftLength[0], lowCornerCell2[1] + _shiftLength[1],
                      lowCornerCell2[2] + _shiftLength[2]};
    const auto [lowCornerCell1, highCornerCell1] = _cellBlocks[_upperLevel]->getCellBoundingBox(offset1);
    // decompose for every level below upperLevel
    for (size_t lowerLevel = 0; lowerLevel < _upperLevel; lowerLevel++) {
      auto startIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(lowCornerCell2);
      auto stopIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(highCornerCell2);
      for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
        for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
          for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
            // Check if the cells are within range of one another
            auto [lowCornerLowerCell, highCornerLowerCell] = _cellBlocks[lowerLevel]->getCellBoundingBox({x, y, z});
            using autopas::utils::ArrayMath::boxDistanceSquared;
            if (boxDistanceSquared(lowCornerLowerCell, highCornerLowerCell, lowCornerCell1, highCornerCell1) >
                _interactionLengthsSquared[lowerLevel]) {
              continue;
            }
            // @todo: Could potentially calculate and use sorting direction in the future
            this->_cellPairOffsetsPerLevel[lowerLevel].emplace_back(
                offset1,
                utils::ThreeDimensionalMapping::threeToOneD(x, y, z,
                                                            _cellBlocks[lowerLevel]->getCellsPerDimensionWithHalo()),
                std::array<double, 3>{0., 0., 0.});
          }
        }
      }
    }
  }

  /**
   * Decomposes the higher-level cell with index cellIndex2 into all lower-level cells inside of it. Then processes
   * the interactions between cell1 and the lower-level cells. The lower-level cells are determined by the bounding
   * box of cell2.
   * @param cell1 The upper-level cell to interact with.
   * @param cellIndex1 The index of the first upper-level cell.
   * @param cellIndex2 The index of the second upper-level cell.
   */
  void decompose2AndProcessCells(ParticleCell_T &cell1, const size_t cellIndex1, const size_t cellIndex2);

  friend class HGC08CellHandlerTest;
};

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void HGC08CellHandler<ParticleCell_T, PairwiseFunctor_T>::decompose2AndProcessCells(ParticleCell_T &cell1,
                                                                                           const size_t cellIndex1,
                                                                                           const size_t cellIndex2) {
  using namespace autopas::utils::ArrayMath::literals;
  // Get bounds of  cell2. Use shifted corners so exact boundary points do not spill into the neighboring cell.
  auto [lowCornerCell2, highCornerCell2] = _cellBlocks[_upperLevel]->getCellBoundingBox(cellIndex2);
  highCornerCell2 = {highCornerCell2[0] - _shiftLength[0], highCornerCell2[1] - _shiftLength[1],
                     highCornerCell2[2] - _shiftLength[2]};

  lowCornerCell2 = {lowCornerCell2[0] + _shiftLength[0], lowCornerCell2[1] + _shiftLength[1],
                    lowCornerCell2[2] + _shiftLength[2]};
  const auto [lowCornerCell1, highCornerCell1] = _cellBlocks[_upperLevel]->getCellBoundingBox(cellIndex1);
  // decompose for every level below upperLevel
  for (size_t lowerLevel = 0; lowerLevel < _upperLevel; lowerLevel++) {
    auto startIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(lowCornerCell2);
    auto stopIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(highCornerCell2);
    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          // Check if the cells are within range of one another
          auto [lowCornerLowerCell, highCornerLowerCell] = _cellBlocks[lowerLevel]->getCellBoundingBox({x, y, z});
          using autopas::utils::ArrayMath::boxDistanceSquared;
          if (boxDistanceSquared(lowCornerLowerCell, highCornerLowerCell, lowCornerCell1, highCornerCell1) >
              _interactionLengthsSquared[lowerLevel]) {
            continue;
          }
          // @todo: Could potentially calculate and use sorting direction in the future
          auto &lowerCell = _cellBlocks[lowerLevel]->getCell({x, y, z});
          this->_cellFunctor.processCellPair(cell1, lowerCell, {0., 0., 0.});
        }
      }
    }
  }
}

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void HGC08CellHandler<ParticleCell_T, PairwiseFunctor_T>::processBaseCell(size_t baseIndexUpperLevel) {
  using namespace autopas::utils::ArrayMath::literals;
  // Intra Level
  for (auto const &[offset1, offset2, r] : this->_cellPairOffsets) {
    const unsigned long cellIndex1 = baseIndexUpperLevel + offset1;
    const unsigned long cellIndex2 = baseIndexUpperLevel + offset2;

    auto &cell1 = _cellBlocks[_upperLevel]->getCell(cellIndex1);
    auto &cell2 = _cellBlocks[_upperLevel]->getCell(cellIndex2);

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
    } else {
      this->_cellFunctor.processCellPair(cell1, cell2, r);
    }
  }
  auto [lowCornerBaseCell, highCornerBaseCell] = _cellBlocks[_upperLevel]->getCellBoundingBox(baseIndexUpperLevel);
  lowCornerBaseCell = lowCornerBaseCell + _shiftLength;

  // Inter Level
  // @todo does it ever make sense to make this parallel/is that possible? Maybe for less thread downtime at the end?
  for (size_t lowerLevel = 0; lowerLevel < _cellPairOffsetsPerLevel.size(); lowerLevel++) {
    const double baseIndexLowerLevel = _cellBlocks[lowerLevel]->get1DIndexOfPosition(lowCornerBaseCell);
    for (auto const &[offset1, offset2, r] : this->_cellPairOffsetsPerLevel[lowerLevel]) {
      const unsigned long cellIndexUpperLevel = baseIndexUpperLevel + offset1;
      const unsigned long cellIndexLowerLevel = baseIndexLowerLevel + offset2;

      auto &cell1 = _cellBlocks[_upperLevel]->getCell(cellIndexUpperLevel);
      auto &cell2 = _cellBlocks[lowerLevel]->getCell(cellIndexLowerLevel);
      this->_cellFunctor.processCellPair(cell1, cell2, r);
    }
  }
}
}  // namespace autopas
