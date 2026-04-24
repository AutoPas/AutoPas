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
 * interactions in the same spatial region.
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
   * @param fittedGrids Whether the grids of the hierarchical grid are fitted to each other.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit HGC08CellHandler(PairwiseFunctor_T *pairwiseFunctor, const std::array<unsigned long, 3> &cellsPerDimension,
                            double interactionLength, const std::array<double, 3> &cellLength,
                            const std::array<unsigned long, 3> &overlap, DataLayoutOption dataLayout, bool useNewton3,
                            const std::vector<internal::CellBlock3D<ParticleCell_T> *> &cellBlocks,
                            const std::vector<double> &interactionLengthsSquared, const size_t upperLevel,
                            bool fittedGrids = true)
      : LCC08CellHandler<ParticleCell_T, PairwiseFunctor_T>(pairwiseFunctor, cellsPerDimension, interactionLength,
                                                            cellLength, overlap, dataLayout, useNewton3),
        _cellBlocks(cellBlocks),
        _interactionLengthsSquared(interactionLengthsSquared),
        _upperLevel(upperLevel),
        _fittedGrids(fittedGrids) {
    using namespace autopas::utils::ArrayMath::literals;
    _shiftLength = _cellBlocks[0]->getCellLength() * 0.01;
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
   * Indicates if the grids are fitted.
   */
  bool _fittedGrids;
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
   * Shift length, which is smaller than the smallest cell width.
   */
  std::array<double, 3> _shiftLength;

  /**
   * Decomposes the higher-level cell with index cellIndex2 into all lower-level cells inside of it. Then processes the
   * interactions between cell1 and the lower-level cells. The lower-level cells are determined by the bounding box of
   * cell2.
   * @param cell1 The upper-level cell to interact with.
   * @param cellIndex1 The index of the first upper-level cell.
   * @param cellIndex2 The index of the second upper-level cell.
   */
  void decompose2AndProcessCells(ParticleCell_T &cell1, const size_t cellIndex1, const size_t cellIndex2);
};

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void HGC08CellHandler<ParticleCell_T, PairwiseFunctor_T>::decompose2AndProcessCells(ParticleCell_T &cell1,
                                                                                           const size_t cellIndex1,
                                                                                           const size_t cellIndex2) {
  using namespace autopas::utils::ArrayMath::literals;

  // get bounds of  cell2
  auto [lowCornerCell2, highCornerCell2] = _cellBlocks[_upperLevel]->getCellBoundingBox(cellIndex2);

  if (_fittedGrids) {
    // Use shifted corners so exact boundary points do not spill into the neighboring cell.
    highCornerCell2 = {highCornerCell2[0] - _shiftLength[0], highCornerCell2[1] - _shiftLength[1],
                       highCornerCell2[2] - _shiftLength[2]};

    lowCornerCell2 = {lowCornerCell2[0] + _shiftLength[0], lowCornerCell2[1] + _shiftLength[1],
                      lowCornerCell2[2] + _shiftLength[2]};
  } else {
    // Replace the upper corner with the lower corner of the diagonally higher corner. Otherwise, we might get ambiguous
    // boundaries, because the high and next low corners don't always give the same result, because of floating point
    // errors.
    auto upperIndex3D = utils::ThreeDimensionalMapping::oneToThreeD(
        cellIndex2, _cellBlocks[_upperLevel]->getCellsPerDimensionWithHalo());
    upperIndex3D = {upperIndex3D[0] + 1, upperIndex3D[1] + 1, upperIndex3D[2] + 1};
    highCornerCell2 = _cellBlocks[_upperLevel]->getCellBoundingBox(upperIndex3D).first;
  }

  const auto [lowCornerCell1, highCornerCell1] = _cellBlocks[_upperLevel]->getCellBoundingBox(cellIndex1);

  // decompose for every level below upperLevel
  for (size_t lowerLevel = 0; lowerLevel < _upperLevel; lowerLevel++) {
    auto startIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(lowCornerCell2);
    auto stopIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(highCornerCell2);

    // if the grids are not fitted, we count lower cells only to the decomposition of an upper cell, if their center
    // lies within, to prevent considering the same pair twice or not at all.
    if (not _fittedGrids) {
      // Calculate the cell centers of the lower-level start and stop index. If they lie outside of the bounds of the
      // upper cell, we shift the start/stop index by one.
      auto [startLow, startHigh] = _cellBlocks[lowerLevel]->getCellBoundingBox(startIndex3D);
      std::array<double, 3> startCellCenter = (startLow + startHigh) * 0.5;
      auto [stopLow, stopHigh] = _cellBlocks[lowerLevel]->getCellBoundingBox(stopIndex3D);
      std::array<double, 3> stopCellCenter = (stopLow + stopHigh) * 0.5;
      if (startCellCenter[0] < lowCornerCell2[0] &&
          startIndex3D[0] < _cellBlocks[lowerLevel]->getCellsPerDimensionWithHalo()[0] - 1) {
        startIndex3D[0]++;
      }
      if (startCellCenter[1] < lowCornerCell2[1] &&
          startIndex3D[1] < _cellBlocks[lowerLevel]->getCellsPerDimensionWithHalo()[1] - 1) {
        startIndex3D[1]++;
      }
      if (startCellCenter[2] < lowCornerCell2[2] &&
          startIndex3D[2] < _cellBlocks[lowerLevel]->getCellsPerDimensionWithHalo()[2] - 1) {
        startIndex3D[2]++;
      }
      if (stopCellCenter[0] >= highCornerCell2[0] && stopIndex3D[0] > 0) {
        stopIndex3D[0]--;
      }
      if (stopCellCenter[1] >= highCornerCell2[1] && stopIndex3D[1] > 0) {
        stopIndex3D[1]--;
      }
      if (stopCellCenter[2] >= highCornerCell2[2] && stopIndex3D[2] > 0) {
        stopIndex3D[2]--;
      }
    }
    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          // Check if the cells are within range of one another
          auto [lowCornerLowerCell, highCornerLowerCell] = _cellBlocks[lowerLevel]->getCellBoundingBox({x, y, z});
          if (lowerLevel < _interactionLengthsSquared.size()) {
            using autopas::utils::ArrayMath::boxDistanceSquared;
            if (boxDistanceSquared(lowCornerLowerCell, highCornerLowerCell, lowCornerCell1, highCornerCell1) >
                _interactionLengthsSquared[lowerLevel]) {
              continue;
            }
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
inline void HGC08CellHandler<ParticleCell_T, PairwiseFunctor_T>::processBaseCell(size_t baseIndex) {
  for (auto const &[offset1, offset2, r] : this->_cellPairOffsets) {
    const unsigned long cellIndex1 = baseIndex + offset1;
    const unsigned long cellIndex2 = baseIndex + offset2;

    auto &cell1 = _cellBlocks[_upperLevel]->getCell(cellIndex1);
    auto &cell2 = _cellBlocks[_upperLevel]->getCell(cellIndex2);
    // if both cells are halo cells, we can skip the interaction,
    if (cell1.getPossibleParticleOwnerships() == OwnershipState::halo and
        cell2.getPossibleParticleOwnerships() == OwnershipState::halo) {
      continue;
    }

    // Intra Level
    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
    } else {
      this->_cellFunctor.processCellPair(cell1, cell2, r);
    }
    // Inter-level only top-down
    if (_upperLevel == 0) {
      continue;
    }

    // Inter Level
    if (not cell1.isEmpty()) {
      decompose2AndProcessCells(cell1, cellIndex1, cellIndex2);
    }
    if (cellIndex1 != cellIndex2 and not cell2.isEmpty()) {
      decompose2AndProcessCells(cell2, cellIndex2, cellIndex1);
    }
  }
}
}  // namespace autopas
