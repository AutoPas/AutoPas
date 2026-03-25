/**
 * @file HGC08CellHandler.h
 * @author S. Seckler
 * @date 10.01.2019
 */

#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

#include "LCC08CellHandler.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
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
class HGC08CellHandler : public LCC08CellHandler<ParticleCell, PairwiseFunctor> {
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
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit HGC08CellHandler(PairwiseFunctor *pairwiseFunctor, const std::array<unsigned long, 3> &cellsPerDimension,
                            double interactionLength, const std::array<double, 3> &cellLength,
                            const std::array<unsigned long, 3> &overlap, DataLayoutOption dataLayout, bool useNewton3,
                            const std::vector<internal::CellBlock3D<ParticleCell> *> &cellBlocks,
                            const std::vector<double> &interactionLengthsSquared, const size_t upperLevel,
                            bool fittedGrids = false)
      : LCC08CellHandler<ParticleCell, PairwiseFunctor>(pairwiseFunctor, cellsPerDimension, interactionLength,
                                                        cellLength, overlap, dataLayout, useNewton3),
        _cellBlocks(cellBlocks),
        _interactionLengthsSquared(interactionLengthsSquared),
        _fittedGrids(fittedGrids),
        _upperLevel(upperLevel) {
    using namespace autopas::utils::ArrayMath::literals;
    _shiftLength = _cellBlocks[0]->getCellLength() * 0.01;
  }

  void processBaseCell(size_t baseIndex);

 protected:
  using Particle = typename ParticleCell::ParticleType;
  using CellBlock = internal::CellBlock3D<ParticleCell>;

  bool _fittedGrids;
  const size_t _upperLevel;
  const std::vector<internal::CellBlock3D<ParticleCell> *> &_cellBlocks;
  const std::vector<double> &_interactionLengthsSquared;
  std::array<double, 3> _shiftLength;
  void decompose2AndProcessCells(ParticleCell &cell1, const size_t cellIndex1, const size_t cellIndex2);

};

template <class ParticleCell, class PairwiseFunctor>
inline void HGC08CellHandler<ParticleCell, PairwiseFunctor>::decompose2AndProcessCells(ParticleCell &cell1,
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
    highCornerCell2 = {std::nextafter(highCornerCell2[0], -std::numeric_limits<double>::infinity()),
                       std::nextafter(highCornerCell2[1], -std::numeric_limits<double>::infinity()),
                       std::nextafter(highCornerCell2[2], -std::numeric_limits<double>::infinity())};
  }

  const auto [lowCell1, highCell1] = _cellBlocks[_upperLevel]->getCellBoundingBox(cellIndex1);
  auto upperIndex3D =
      utils::ThreeDimensionalMapping::oneToThreeD(cellIndex2, _cellBlocks[_upperLevel]->getCellsPerDimensionWithHalo());

  // decompose for every level below upperLevel
  for (size_t lowerLevel = 0; lowerLevel < _upperLevel; lowerLevel++) {
    // see if all of this is worth it
    auto startIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(lowCornerCell2);
    auto stopIndex3D = _cellBlocks[lowerLevel]->get3DIndexOfPosition(highCornerCell2);

    // if the grids are not fitted, we count lower cells only to the decomposition of an upper cell, if their center
    // lies within, to prevent considering the same pair twice
    if (!_fittedGrids) {
      const auto [low, high] = _cellBlocks[lowerLevel]->getCellBoundingBox(startIndex3D);
      std::array<double, 3> startCellCenter = {0.5 * (low[0] + high[0]), 0.5 * (low[1] + high[1]),
                                               0.5 * (low[2] + high[2])};
      const auto [stopLow, stopHigh] = _cellBlocks[lowerLevel]->getCellBoundingBox(stopIndex3D);
      std::array<double, 3> stopCellCenter = {0.5 * (stopLow[0] + stopHigh[0]), 0.5 * (stopLow[1] + stopHigh[1]),
                                              0.5 * (stopLow[2] + stopHigh[2])};
      startCellCenter = {0.5 * (lowCornerCell2[0] + highCornerCell2[0]), 0.5 * (lowCornerCell2[1] + highCornerCell2[1]),
                         0.5 * (lowCornerCell2[2] + highCornerCell2[2])};
      // do pretty conservative checks here, so we dont miss anything
      // double checks are prevented by ownership checks in the loop
      if (startCellCenter[0] + _shiftLength[0] < lowCornerCell2[0] &&
          startIndex3D[0] < _cellBlocks[lowerLevel]->getCellsPerDimensionWithHalo()[0] - 1) {
        startIndex3D[0]++;
      }
      if (startCellCenter[1] + _shiftLength[1] < lowCornerCell2[1] &&
          startIndex3D[1] < _cellBlocks[lowerLevel]->getCellsPerDimensionWithHalo()[1] - 1) {
        startIndex3D[1]++;
      }
      if (startCellCenter[2] + _shiftLength[2] < lowCornerCell2[2] &&
          startIndex3D[2] < _cellBlocks[lowerLevel]->getCellsPerDimensionWithHalo()[2] - 1) {
        startIndex3D[2]++;
      }
      // I think 0 checks are not necessary, but eh
      if (stopCellCenter[0] - _shiftLength[0] > highCornerCell2[0] && stopIndex3D[0] > 0) {
        stopIndex3D[0]--;
      }
      if (stopCellCenter[1] - _shiftLength[1] > highCornerCell2[1] && stopIndex3D[1] > 0) {
        stopIndex3D[1]--;
      }
      if (stopCellCenter[2] - _shiftLength[2] > highCornerCell2[2] && stopIndex3D[2] > 0) {
        stopIndex3D[2]--;
      }
    }
    // debug info start
    // const auto startCellBoundingBoxSecond = _cellBlocks[lowerLevel]->getCellBoundingBox(startIndex3D);
    // const auto stopCellBoundingBoxSecond = _cellBlocks[lowerLevel]->getCellBoundingBox(stopIndex3D);
    /*    std::stringstream debugStream;
    debugStream << "cell2 low: [" << lowCornerCell2[0] << ", " << lowCornerCell2[1] << ", " << lowCornerCell2[2]
                << "] high: [" << highCornerCell2[0] << ", " << highCornerCell2[1] << ", " << highCornerCell2[2]
                << "] |\n "
                << "lowerLevel " << lowerLevel << " second-start low: [" << startCellBoundingBoxSecond.first[0] << ", "
                << startCellBoundingBoxSecond.first[1] << ", " << startCellBoundingBoxSecond.first[2] << "] high: ["
                << startCellBoundingBoxSecond.second[0] << ", " << startCellBoundingBoxSecond.second[1] << ", "
                << startCellBoundingBoxSecond.second[2]
                << "] |\n "
                   " second-stop low: ["
                << stopCellBoundingBoxSecond.first[0] << ", " << stopCellBoundingBoxSecond.first[1] << ", "
                << stopCellBoundingBoxSecond.first[2] << "] high: [" << stopCellBoundingBoxSecond.second[0] << ", "
                << stopCellBoundingBoxSecond.second[1] << ", " << stopCellBoundingBoxSecond.second[2] << "]\n";

    std::cout << debugStream.str() << std::endl;
*/
    // debug info end

    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          using autopas::utils::ArrayMath::boxDistanceSquared;
          // check if number is in range
          /*          auto [low, high] = _cellBlocks[lowerLevel]->getCellBoundingBox({x, y, z});
          if (lowerLevel < _interactionLengthsSquared.size()) {
            if (boxDistanceSquared(low, high, lowCell1, highCell1) > _interactionLengthsSquared[lowerLevel]) {
              continue;
            }
          }
*/
          if (!_fittedGrids && (x == stopIndex3D[0] || y == stopIndex3D[1] || z == stopIndex3D[2] ||
                                x == startIndex3D[0] || y == startIndex3D[1] || z == startIndex3D[2])) {
            // Canonical ownership for non-fitted grids:
            // assign lower cell only to the upper cell that contains its center.
            std::array<double, 3> lowerCenter{0.5 * (low[0] + high[0]), 0.5 * (low[1] + high[1]),
                                              0.5 * (low[2] + high[2])};
            lowerCenter = {std::nextafter(lowerCenter[0], std::numeric_limits<double>::infinity()),
                           std::nextafter(lowerCenter[1], std::numeric_limits<double>::infinity()),
                           std::nextafter(lowerCenter[2], std::numeric_limits<double>::infinity())};
            auto ownerUpperIndex3D = _cellBlocks[_upperLevel]->get3DIndexOfPosition(lowerCenter);
            if (ownerUpperIndex3D != upperIndex3D) {
              continue;
            }
          }

          // Could potentially calculate and use sorting direction in the future
          auto &lowerCell = _cellBlocks[lowerLevel]->getCell({x, y, z});
          this->_cellFunctor.processCellPair(cell1, lowerCell, {0., 0., 0.});
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor>
inline void HGC08CellHandler<ParticleCell, PairwiseFunctor>::processBaseCell(size_t baseIndex) {
  for (auto const &[offset1, offset2, r] : this->_cellPairOffsets) {
    const unsigned long cellIndex1 = baseIndex + offset1;
    const unsigned long cellIndex2 = baseIndex + offset2;

    auto &cell1 = _cellBlocks[_upperLevel]->getCell(cellIndex1);
    auto &cell2 = _cellBlocks[_upperLevel]->getCell(cellIndex2);
    // if both cells are halo cells, we can skip the interaction,
    if (cell1.getPossibleParticleOwnerships() == OwnershipState::halo &&
        cell2.getPossibleParticleOwnerships() == OwnershipState::halo) {
      continue;
    }

    // Intra Level
    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
      decompose2AndProcessCells(cell1, cellIndex1, cellIndex2);
    } else {
      // Intra Level
      this->_cellFunctor.processCellPair(cell1, cell2, r);
      if (_upperLevel == 0) {
        continue;
      }

      // Inter Level
      if (!cell1.isEmpty()) {
        decompose2AndProcessCells(cell1, cellIndex1, cellIndex2);
      }
      if (!cell2.isEmpty()) {
        decompose2AndProcessCells(cell2, cellIndex2, cellIndex1);
      }
    }
  }
}
}  // namespace autopas
