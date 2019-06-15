/**
 * @file C04CellHandler.h
 * @author S. Seckler
 * @date 10.01.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c04 base step.
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
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C04CellHandler {
 public:
  /**
   * Constructor of the c04 traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   */
  explicit C04CellHandler(PairwiseFunctor *pairwiseFunctor, std::array<unsigned long, 3> cellsPerDimension,
                          const double cutoff = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0},
                          const std::array<unsigned long, 3> &overlap = {1ul, 1ul, 1ul})
      : _cellFunctor(pairwiseFunctor),
        _cellPairOffsets{},
        _cutoff(cutoff),
        _cellLength(cellLength),
        _overlap(overlap),
        _cellsPerDimension(cellsPerDimension) {
    computeOffsets(cellsPerDimension);
  }

  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner (=base index) of a 2x2x2 block of cells.
   * @param cells vector of all cells.
   * @param baseIndex Index respective to which box is constructed.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  /**
   * Resize all buffers to match the current number of threads.
   */
  void resizeBuffers();

 private:
  /**
   * Computes pairs for the block used in processBaseCell().
   * The algorithm used to generate the cell pairs can be visualized with a python script, which can be found in
   * docs/C08TraversalScheme.py
   * @param cellsPerDimension
   */
  void computeOffsets(std::array<unsigned long, 3> cellsPerDimension);

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, useNewton3> _cellFunctor;

  /**
   * Pair sets for processBaseCell().
   */
  std::vector<std::vector<std::pair<unsigned long, unsigned long>>> _cellPairOffsets;

  /**
   * cutoff radius.
   */
  const double _cutoff;

  /**
   * cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;

  /**
   * overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  const std::array<unsigned long, 3> _overlap;

  std::array<unsigned long, 3> _cellsPerDimension;

  std::vector<std::vector<int>> baseOffsets;

  std::vector<std::vector<ParticleCell>> _combinationSlices;
  std::vector<unsigned int> _currentSlices;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04CellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  using std::pair;

  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, _cellsPerDimension);

  const auto threadID = autopas_get_thread_num();
  auto &currentSlice = _currentSlices[threadID];
  auto &combinationSlice = _combinationSlices[threadID];

  // First cell needs to initialize whole buffer
  if (x == this->_overlap[0]) {
    currentSlice = 0;
    for (unsigned int offsetSlice = 0; offsetSlice < baseOffsets.size(); offsetSlice++) {
      combinationSlice[offsetSlice]._particleSoABuffer.clear();
      for (const auto offset : baseOffsets[offsetSlice]) {
        const unsigned long otherIndex = baseIndex + offset;
        ParticleCell &otherCell = cells[otherIndex];
        combinationSlice[offsetSlice]._particleSoABuffer.append(otherCell._particleSoABuffer);
      }
    }
  } else {
    combinationSlice[currentSlice]._particleSoABuffer.clear();

    for (const auto offset : baseOffsets.back()) {
      const unsigned long otherIndex = baseIndex + offset;
      ParticleCell &otherCell = cells[otherIndex];
      combinationSlice[currentSlice]._particleSoABuffer.append(otherCell._particleSoABuffer);
    }

    ++currentSlice %= baseOffsets.size();
  }

  for (auto &slice : _cellPairOffsets) {
    for (auto &current_pair : slice) {
      unsigned long offset1 = current_pair.first;
      unsigned long cellIndex1 = baseIndex + offset1;

      unsigned long offset2 = current_pair.second;
      unsigned long cellIndex2 = baseIndex + offset2;

      ParticleCell &cell1 = cells[cellIndex1];
      ParticleCell &cell2 = cells[cellIndex2];

      if (cellIndex1 == cellIndex2) {
        this->_cellFunctor.processCell(cell1);
      } else {
        this->_cellFunctor.processCellPair(cell1, cell2);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04CellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::computeOffsets(
    std::array<unsigned long, 3> cellsPerDimension) {
  using std::make_pair;

  //////////////////////////////
  // @TODO: Replace following lines with vector to support asymmetric cells
  const unsigned long ov1 = _overlap[0] + 1;
  const unsigned long ov1_squared = ov1 * ov1;
  //////////////////////////////

  baseOffsets.resize(ov1);

  std::array<unsigned long, 3> overlap_1 = ArrayMath::addScalar(_overlap, 1ul);

  std::vector<unsigned long> cellOffsets;
  cellOffsets.reserve(overlap_1[0] * overlap_1[1] * overlap_1[2]);

  _cellPairOffsets.clear();
  _cellPairOffsets.resize(overlap_1[0]);

  const auto cutoffSquare(this->_cutoff * this->_cutoff);

  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
      for (unsigned long z = 0ul; z <= _overlap[2]; ++z) {
        cellOffsets.push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension));
      }
    }
  }
  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
      baseOffsets[x].push_back(ov1_squared * x + ov1 * y);
      for (unsigned long z = 0ul; z <= _overlap[2]; ++z) {
        const unsigned long offset = cellOffsets[ov1_squared * x + ov1 * y];
        // origin
        {
          // check whether cell is within cutoff radius
          auto distVec =
              ArrayMath::mul({std::max(0.0, x - 1.0), std::max(0.0, y - 1.0), std::max(0.0, z - 1.0)}, _cellLength);
          const auto distSquare = ArrayMath::dot(distVec, distVec);
          if (distSquare <= cutoffSquare) {
            _cellPairOffsets[x].push_back(make_pair(cellOffsets[z], offset));
          }
        }
        // back left
        if (y != _overlap[1] and z != 0) {
          // check whether cell is within cutoff radius
          auto distVec = ArrayMath::mul(
              {std::max(0.0, x - 1.0), std::max(0.0, _overlap[1] - y - 1.0), std::max(0.0, z - 1.0)}, _cellLength);
          const auto distSquare = ArrayMath::dot(distVec, distVec);
          if (distSquare <= cutoffSquare) {
            _cellPairOffsets[x].push_back(make_pair(cellOffsets[ov1_squared - ov1 + z], offset));
          }
        }
        // front right
        if (x != _overlap[0] and (y != 0 or z != 0)) {
          // check whether cell is within cutoff radius
          auto distVec = ArrayMath::mul(
              {std::max(0.0, _overlap[0] - x - 1.0), std::max(0.0, y - 1.0), std::max(0.0, z - 1.0)}, _cellLength);
          const auto distSquare = ArrayMath::dot(distVec, distVec);
          if (distSquare <= cutoffSquare) {
            _cellPairOffsets[x].push_back(make_pair(cellOffsets[ov1_squared * _overlap[0] + z], offset));
          }
        }
        // back right
        if (y != _overlap[1] and x != _overlap[0] and z != 0) {
          // check whether cell is within cutoff radius
          auto distVec = ArrayMath::mul(
              {std::max(0.0, _overlap[0] - x - 1.0), std::max(0.0, _overlap[1] - y - 1.0), std::max(0.0, z - 1.0)},
              _cellLength);
          const auto distSquare = ArrayMath::dot(distVec, distVec);
          if (distSquare <= cutoffSquare) {
            _cellPairOffsets[x].push_back(make_pair(cellOffsets[ov1_squared * ov1 - ov1 + z], offset));
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void C04CellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::resizeBuffers() {
  const unsigned int numThreads = static_cast<unsigned int>(autopas_get_max_threads());
  if (_combinationSlices.size() != numThreads) {
    _combinationSlices.resize(numThreads);
    const auto cellOffsetsSize = _overlap[0] + 1;
    std::for_each(_combinationSlices.begin(), _combinationSlices.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _currentSlices.resize(numThreads);
  }
}

}  // namespace autopas
