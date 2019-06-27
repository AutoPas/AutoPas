/**
 * @file C04SoACellHandler.h
 * @author C.Menges
 * @date 04.06.2019
 */

#pragma once

#include <error.h>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step, but rather use 4 instead of 8 colors for domain
 * coloring. Using 4 instead of 8 colors allows combination of SoA buffers.
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
class C04SoACellHandler {
 public:
  /**
   * Constructor of the c04 traversal with combined SoA buffers.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   */
  explicit C04SoACellHandler(PairwiseFunctor *pairwiseFunctor, std::array<unsigned long, 3> cellsPerDimension,
                             const double cutoff = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0},
                             const std::array<unsigned long, 3> &overlap = {1ul, 1ul, 1ul})
      : _cellFunctor(pairwiseFunctor),
        _cutoff(cutoff),
        _cellLength(cellLength),
        _overlap(overlap),
        _cellsPerDimension(cellsPerDimension),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(unsigned int)) {
    computeOffsets(cellsPerDimension);
  }

  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner of a 2x2x2 block of cells.
   * @param cells vector of all cells.
   * @param x cell index x-axis
   * @param y cell index y-axis
   * @param z cell index z-axis
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
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, useNewton3>
      _cellFunctor;

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

  const std::array<unsigned long, 3> _cellsPerDimension;

  /**
   * Cells containing combined SoA buffers
   */
  std::vector<std::vector<ParticleCell>> _combinationSlices;

  /**
   * Current position in circular buffer for each thread
   */
  std::vector<unsigned int> _currentSlices;

  /**
   * cell offsets of "base plate"
   */
  std::vector<std::vector<unsigned long>> _baseOffsets;

  /**
   * Type to store an interval of indices.
   */
  using interval_t = std::pair<unsigned long, unsigned long>;

  /**
   * Interactions cell <-> baseplate intervall
   */
  std::vector<std::vector<std::pair<unsigned long, interval_t>>> _offsets;

  /**
   * Partial sums of sizes of combined buffers to determine start and end quickly
   */
  std::vector<std::vector<std::vector<unsigned long>>> _combinationSlicesOffsets;

  /**
   * Offset factor to avoid false sharing.
   */
  const unsigned int _cacheOffset;

  /**
   * Writes buffer content back to cell.
   * @param cells
   * @param baseIndex Index of base cell.
   * @param combinationSlice
   * @param combinationSlicesOffsets
   * @param bufferSlice Index of slice in combinationSlice (source)
   * @param cellSlice Index of slice in _baseOffsets (destination)
   */
  void writeBufferIntoCell(std::vector<ParticleCell> &cells, const unsigned long baseIndex,
                           std::vector<ParticleCell> &combinationSlice,
                           std::vector<std::vector<unsigned long>> &combinationSlicesOffsets,
                           const unsigned long bufferSlice, const unsigned long cellSlice);

  /**
   * Writes cell content into buffer.
   * @param cells
   * @param baseIndex Index of base cell.
   * @param combinationSlice
   * @param combinationSlicesOffsets
   * @param bufferSlice Index of slice in combinationSlice (destination)
   * @param cellSlice Index of slice in _baseOffsets (source)
   */
  void writeCellIntoBuffer(const std::vector<ParticleCell> &cells, const unsigned long baseIndex,
                           std::vector<ParticleCell> &combinationSlice,
                           std::vector<std::vector<unsigned long>> &combinationSlicesOffsets,
                           const unsigned int bufferSlice, const unsigned int cellSlice);

  /**
   * Creates offset intervals (stored in offset) from _cellPairOffsets
   * @param _cellPairOffsets Source for interval creation.
   */
  void setupIntervals(std::vector<std::vector<std::pair<unsigned long, unsigned long>>> &_cellPairOffsets);
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04SoACellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, _cellsPerDimension);

  // get all information for current thread
  const auto threadID = autopas_get_thread_num();
  auto &currentSlice = _currentSlices[threadID * _cacheOffset];
  auto &combinationSlice = _combinationSlices[threadID];
  auto &combinationSlicesOffsets = _combinationSlicesOffsets[threadID];

  const size_t numSlices = _baseOffsets.size();

  // First cell needs to initialize whole buffer
  if (x == 0ul) {
    currentSlice = 0ul;
    for (unsigned int slice = 0ul; slice < numSlices; slice++) {
      writeCellIntoBuffer(cells, baseIndex, combinationSlice, combinationSlicesOffsets, slice, slice);
    }
  } else {
    writeCellIntoBuffer(cells, baseIndex, combinationSlice, combinationSlicesOffsets, currentSlice, numSlices - 1ul);

    ++currentSlice %= numSlices;
  }

  // compute interactions
  for (unsigned long slice = 0; slice < numSlices; slice++) {
    for (const auto &current_pair : _offsets[(slice + currentSlice) % numSlices]) {
      const unsigned long offset1 = current_pair.first;

      ParticleCell *cell1;

      const auto interval = current_pair.second;
      // special cases (cell1 one is also stored in a combination slice)
      // base cell
      if (offset1 == 0ul) {
        // two cases intervall in current  stripe
        cell1 = &combinationSlice[currentSlice];
        cell1->_particleSoABuffer.setViewStart(0);
        if (slice == currentSlice) {
          // process stripe with itself
          // make sure no previously applied view is active
          cell1->_particleSoABuffer.setViewLength(-1);
          this->_cellFunctor.processCell(*cell1);
          continue;
        } else {
          // interval in other stripe
          cell1->_particleSoABuffer.setViewLength(cells[baseIndex].numParticles());
        }
      } else if (offset1 == _baseOffsets.front().back()) {
        cell1 = &combinationSlice[currentSlice];
        cell1->_particleSoABuffer.setViewStart(
            combinationSlicesOffsets[slice][combinationSlicesOffsets[slice].size() - 2]);
        cell1->_particleSoABuffer.setViewLength(-1l);
      } else {
        const unsigned long cellIndex1 = baseIndex + offset1;
        cell1 = &cells[cellIndex1];
      }

      auto &currentCS = combinationSlice[slice];
      const auto &currentCSOffsets = combinationSlicesOffsets[slice];
      // set view start
      currentCS._particleSoABuffer.setViewStart(currentCSOffsets[interval.first]);
      // set view length
      const auto length = currentCSOffsets[interval.second] - currentCSOffsets[interval.first];
      currentCS._particleSoABuffer.setViewLength(length);

      // process cell pair
      this->_cellFunctor.processCellPair(*cell1, currentCS);
    }
  }

  // write information of combined SoA buffers back to cell
  writeBufferIntoCell(cells, baseIndex, combinationSlice, combinationSlicesOffsets, currentSlice, 0ul);

  // last cell in stripe need to save whole buffer
  if (x == _cellsPerDimension[0] - _overlap[0] - 1l) {
    for (unsigned long slice = 1; slice < numSlices; slice++) {
      const long bufferSlice = (currentSlice + slice) % numSlices;
      writeBufferIntoCell(cells, baseIndex, combinationSlice, combinationSlicesOffsets, bufferSlice, slice);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04SoACellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::writeCellIntoBuffer(
    const std::vector<ParticleCell> &cells, const unsigned long baseIndex, std::vector<ParticleCell> &combinationSlice,
    std::vector<std::vector<unsigned long>> &combinationSlicesOffsets, const unsigned int bufferSlice,
    const unsigned int cellSlice) {
  // delete old data
  combinationSlice[bufferSlice]._particleSoABuffer.clear();
  combinationSlicesOffsets[bufferSlice].clear();
  // fill buffer and compute partial sums
  unsigned long sum = 0ul;
  combinationSlicesOffsets[bufferSlice].push_back(sum);
  for (const auto offset : _baseOffsets[cellSlice]) {
    const unsigned long otherIndex = baseIndex + offset;
    const ParticleCell &otherCell = cells[otherIndex];
    combinationSlice[bufferSlice]._particleSoABuffer.append(otherCell._particleSoABuffer);
    sum += otherCell.numParticles();
    combinationSlicesOffsets[bufferSlice].push_back(sum);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04SoACellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::writeBufferIntoCell(
    std::vector<ParticleCell> &cells, const unsigned long baseIndex, std::vector<ParticleCell> &combinationSlice,
    std::vector<std::vector<unsigned long>> &combinationSlicesOffsets, const unsigned long bufferSlice,
    const unsigned long cellSlice) {
  auto &buffer = combinationSlice[bufferSlice]._particleSoABuffer;
  // Make sure no view length is set
  buffer.setViewLength(-1l);
  for (long i = _baseOffsets[cellSlice].size() - 1l; i >= 0; i--) {
    const auto start = combinationSlicesOffsets[bufferSlice][i];
    const auto end = combinationSlicesOffsets[bufferSlice][i + 1];

    if (start == end) {
      continue;
    }
    buffer.setViewStart(0);
    buffer.resizeArrays(end);
    buffer.setViewStart(start);

    const unsigned long currentOffset = baseIndex + _baseOffsets[cellSlice][i];
    // clear old cell buffer
    cells[currentOffset]._particleSoABuffer.clear();
    // make sure everything is correct
    if (buffer.getNumParticles() != cells[currentOffset].numParticles()) {
      const auto pos = utils::ThreeDimensionalMapping::oneToThreeD(currentOffset, _cellsPerDimension);
      AutoPasLog(error,
                 "Particle number in SoA buffer and cell doesn't match. current position: [{} {} {}] is: {} should: {}",
                 pos[0], pos[1], pos[2], buffer.getNumParticles(), cells[currentOffset].numParticles());
    }
    // append new cell buffer
    cells[currentOffset]._particleSoABuffer.append(buffer);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04SoACellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::computeOffsets(
    std::array<unsigned long, 3> cellsPerDimension) {
  using std::make_pair;

  std::vector<std::vector<std::pair<unsigned long, unsigned long>>> _cellPairOffsets;

  //////////////////////////////
  // @TODO: Replace following lines with vector to support asymmetric cells
  const unsigned long ov1 = _overlap[0] + 1;
  const unsigned long ov1_squared = ov1 * ov1;
  //////////////////////////////

  _baseOffsets.resize(ov1);

  std::array<unsigned long, 3> overlap_1 = ArrayMath::addScalar(_overlap, 1ul);

  std::vector<unsigned long> cellOffsets;
  cellOffsets.reserve(overlap_1[0] * overlap_1[1] * overlap_1[2]);

  _cellPairOffsets.clear();
  _cellPairOffsets.resize(overlap_1[0]);

  const auto cutoffSquare(this->_cutoff * this->_cutoff);

  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
      _baseOffsets[x].push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, 0ul, cellsPerDimension));
      for (unsigned long z = 0ul; z <= _overlap[2]; ++z) {
        cellOffsets.push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension));
      }
    }
  }
  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
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
  setupIntervals(_cellPairOffsets);
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04SoACellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::setupIntervals(
    std::vector<std::vector<std::pair<unsigned long, unsigned long>>> &_cellPairOffsets) {
  // Create intervals
  const unsigned long numStripes = _cellPairOffsets.size();
  _offsets.resize(numStripes);

  // iterate over all stripes
  for (unsigned long i = 0; i < numStripes; ++i) {
    auto &currentStripe = _cellPairOffsets[i];
    // Sort
    std::stable_sort(currentStripe.begin(), currentStripe.end(),
                     [](const auto &a, const auto &b) -> bool { return a.first < b.first; });

    // Collect intervals
    unsigned long current = currentStripe.front().first;
    unsigned long startID =
        std::distance(_baseOffsets[i].begin(),
                      std::find(_baseOffsets[i].begin(), _baseOffsets[i].end(), currentStripe.front().second));
    unsigned long endID = startID;
    for (unsigned long j = 0; j < currentStripe.size(); j++) {
      if (current != currentStripe[j].first) {
        auto interval = std::make_pair(startID, endID);
        _offsets[i].push_back(std::make_pair(current, interval));
        startID = std::distance(_baseOffsets[i].begin(),
                                std::find(_baseOffsets[i].begin(), _baseOffsets[i].end(), currentStripe[j].second));
        endID = startID;
        current = currentStripe[j].first;
      }
      endID++;
    }
    // last interval
    auto interval = std::make_pair(startID, endID);
    _offsets[i].push_back(std::make_pair(current, interval));
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void C04SoACellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::resizeBuffers() {
  const unsigned int numThreads = static_cast<unsigned int>(autopas_get_max_threads());
  if (_combinationSlices.size() != numThreads) {
    _combinationSlices.resize(numThreads);
    const auto cellOffsetsSize = _overlap[0] + 1;
    std::for_each(_combinationSlices.begin(), _combinationSlices.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _combinationSlicesOffsets.resize(numThreads);
    std::for_each(_combinationSlicesOffsets.begin(), _combinationSlicesOffsets.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _currentSlices.resize(numThreads * _cacheOffset);
  }
}

}  // namespace autopas
