/**
 * @file LCC08CellHandler.h
 * @author S. Seckler
 * @date 10.01.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

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
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC08CellHandler {
 public:
  /**
   * Constructor of the LCC08CellHandler.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC08CellHandler(PairwiseFunctor *pairwiseFunctor, std::array<uint64_t, 3> cellsPerDimension,
                            const double interactionLength, const std::array<double, 3> &cellLength,
                            const std::array<uint64_t, 3> &overlap)
      : _cellFunctor(pairwiseFunctor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/),
        _cellPairOffsets{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlap(overlap) {
    computeOffsets(cellsPerDimension);
  }

  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner (=base index) of a 2x2x2 block of cells.
   * @param cells vector of all cells.
   * @param baseIndex Index respective to which box is constructed.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, uint64_t baseIndex);

 private:
  /**
   * Computes pairs for the block used in processBaseCell().
   * The algorithm used to generate the cell pairs can be visualized with a python script, which can be found in
   * docs/C08TraversalScheme.py
   * @param cellsPerDimension
   */
  void computeOffsets(std::array<uint64_t, 3> cellsPerDimension);

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, dataLayout, useNewton3>
      _cellFunctor;

  /**
   * Pair sets for processBaseCell().
   * Values are: offset of first cell, offset of second cell, sorting direction.
   */
  std::vector<std::tuple<uint64_t, uint64_t, std::array<double, 3>>> _cellPairOffsets;

  /**
   * Interaction length (cutoff + skin).
   */
  const double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;

  /**
   * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  const std::array<uint64_t, 3> _overlap;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, uint64_t baseIndex) {
  for (auto const &[offset1, offset2, r] : _cellPairOffsets) {
    const uint64_t cellIndex1 = baseIndex + offset1;
    const uint64_t cellIndex2 = baseIndex + offset2;

    ParticleCell &cell1 = cells[cellIndex1];
    ParticleCell &cell2 = cells[cellIndex2];

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
    } else {
      this->_cellFunctor.processCellPair(cell1, cell2, r);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::computeOffsets(
    std::array<uint64_t, 3> cellsPerDimension) {
  using std::make_pair;

  //////////////////////////////
  // @TODO: Replace following lines with vector to support asymmetric cells
  const uint64_t ov1 = _overlap[0] + 1;
  const uint64_t ov1_squared = ov1 * ov1;
  //////////////////////////////

  std::array<uint64_t, 3> overlap_1 = utils::ArrayMath::addScalar(_overlap, 1ull);

  std::vector<uint64_t> cellOffsets;
  cellOffsets.reserve(overlap_1[0] * overlap_1[1] * overlap_1[2]);

  _cellPairOffsets.clear();

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  // constants 0.0 and 1.0 with correct type
  const double zero = 0.0;
  const double one = 1.0;

  for (uint64_t x = 0ul; x <= _overlap[0]; ++x) {
    for (uint64_t y = 0ul; y <= _overlap[1]; ++y) {
      for (uint64_t z = 0ul; z <= _overlap[2]; ++z) {
        cellOffsets.push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension));
      }
    }
  }
  for (uint64_t x = 0ul; x <= _overlap[0]; ++x) {
    for (uint64_t y = 0ul; y <= _overlap[1]; ++y) {
      for (uint64_t z = 0ul; z <= _overlap[2]; ++z) {
        const uint64_t offset = cellOffsets[ov1_squared * x + ov1 * y];
        // origin
        {
          // check whether cell is within interaction length
          auto distVec = utils::ArrayMath::mul(
              {std::max(zero, x - one), std::max(zero, y - one), std::max(zero, z - one)}, _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            _cellPairOffsets.push_back(std::make_tuple(cellOffsets[z], offset, utils::ArrayMath::normalize(distVec)));
          }
        }
        // back left
        if (y != _overlap[1] and z != 0) {
          // check whether cell is within interaction length
          auto distVec = utils::ArrayMath::mul(
              {std::max(zero, x - one), std::max(zero, _overlap[1] - y - one), std::max(zero, z - one)}, _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            _cellPairOffsets.push_back(
                std::make_tuple(cellOffsets[ov1_squared - ov1 + z], offset, utils::ArrayMath::normalize(distVec)));
          }
        }
        // front right
        if (x != _overlap[0] and (y != 0 or z != 0)) {
          // check whether cell is within interaction length
          auto distVec = utils::ArrayMath::mul(
              {std::max(zero, _overlap[0] - x - one), std::max(zero, y - one), std::max(zero, z - one)}, _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            _cellPairOffsets.push_back(std::make_tuple(cellOffsets[ov1_squared * _overlap[0] + z], offset,
                                                       utils::ArrayMath::normalize(distVec)));
          }
        }
        // back right
        if (y != _overlap[1] and x != _overlap[0] and z != 0) {
          // check whether cell is within interaction length
          auto distVec = utils::ArrayMath::mul(
              {std::max(zero, _overlap[0] - x - one), std::max(zero, _overlap[1] - y - one), std::max(zero, z - one)},
              _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            _cellPairOffsets.push_back(std::make_tuple(cellOffsets[ov1_squared * ov1 - ov1 + z], offset,
                                                       utils::ArrayMath::normalize(distVec)));
          }
        }
      }
    }
  }
}

}  // namespace autopas
