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
        _cutoff(cutoff),
        _cellLength(cellLength),
        _overlap(overlap),
        _cellsPerDimension(cellsPerDimension) {
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

  std::array<unsigned long, 3> _cellsPerDimension;

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
   * Interactions cell <-> baseplate intervall
   */
  std::vector<std::vector<std::pair<unsigned long, std::pair<unsigned long, unsigned long>>>> _offsets;

  /**
   * Partial sums of sizes of combined buffers to determine start and end quickly
   */
  std::vector<std::vector<std::vector<unsigned long>>> _combinationSlicesOffsets;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04CellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, _cellsPerDimension);

  // get all information for current thread
  const auto threadID = autopas_get_thread_num();
  auto &currentSlice = _currentSlices[threadID];
  auto &combinationSlice = _combinationSlices[threadID];
  auto &combinationSlicesOffsets = _combinationSlicesOffsets[threadID];

  // First cell needs to initialize whole buffer
  if (x == 0ul) {
    currentSlice = 0;
    for (unsigned int offsetSlice = 0; offsetSlice < _baseOffsets.size(); offsetSlice++) {
      // delete old data
      combinationSlice[offsetSlice]._particleSoABuffer.clear();
      combinationSlicesOffsets[offsetSlice].clear();
      // fill buffer and compute partial sums
      unsigned long sum = 0ul;
      combinationSlicesOffsets[offsetSlice].push_back(sum);
      for (const auto offset : _baseOffsets[offsetSlice]) {
        const unsigned long otherIndex = baseIndex + offset;
        ParticleCell &otherCell = cells[otherIndex];
        combinationSlice[offsetSlice]._particleSoABuffer.append(otherCell._particleSoABuffer);
        sum += otherCell.numParticles();
        combinationSlicesOffsets[offsetSlice].push_back(sum);
      }
    }
  } else {
    // delete old data
    combinationSlice[currentSlice]._particleSoABuffer.clear();
    combinationSlicesOffsets[currentSlice].clear();
    // fill first stripe buffer and compute partial sums
    unsigned long sum = 0ul;
    combinationSlicesOffsets[currentSlice].push_back(sum);
    for (const auto offset : _baseOffsets.back()) {
      const unsigned long otherIndex = baseIndex + offset;
      ParticleCell &otherCell = cells[otherIndex];
      combinationSlice[currentSlice]._particleSoABuffer.append(otherCell._particleSoABuffer);
      sum += otherCell.numParticles();
      combinationSlicesOffsets[currentSlice].push_back(sum);
    }

    ++currentSlice %= _baseOffsets.size();
  }

  // compute interactions
  for (unsigned long slice = 0; slice < combinationSlicesOffsets.size(); slice++) {
    for (auto &current_pair : _offsets[slice]) {
      unsigned long offset1 = current_pair.first;
      unsigned long cellIndex1 = baseIndex + offset1;

      ParticleCell *cell1 = &cells[cellIndex1];

      auto interval = current_pair.second;
      // special cases front left/backleft
      // base cell
      if (offset1 == _baseOffsets.front().front()) {
        // two cases intervall in current  stripe
        cell1 = &combinationSlice[currentSlice % combinationSlice.size()];
        if (slice == 0ul) {
          // process stripe with itself
          // make sure no previously applied view is active
          cell1->_particleSoABuffer.setViewStart(0);
          cell1->_particleSoABuffer.setViewLength(-1);
          this->_cellFunctor.processCell(*cell1);
          continue;
        } else {
          // interval in other stripe
          cell1->_particleSoABuffer.setViewStart(0);
          cell1->_particleSoABuffer.setViewLength(cells[cellIndex1].numParticles());
        }
      } else if (offset1 == _baseOffsets.front().back()) {
        cell1 = &combinationSlice[currentSlice % combinationSlice.size()];
        cell1->_particleSoABuffer.setViewStart(
            combinationSlicesOffsets[slice][combinationSlicesOffsets[slice].size() - 2]);
        cell1->_particleSoABuffer.setViewLength(-1l);
      }

      auto &currentCombinationSlice = combinationSlice[(currentSlice + slice + 1) % combinationSlice.size()];
      // set view start
      currentCombinationSlice._particleSoABuffer.setViewStart(combinationSlicesOffsets[slice][interval.first]);
      // set view length
      currentCombinationSlice._particleSoABuffer.setViewLength(combinationSlicesOffsets[slice][interval.second] -
                                                               combinationSlicesOffsets[slice][interval.first]);

      // process cell pair
      this->_cellFunctor.processCellPair(*cell1, currentCombinationSlice);
    }
  }

  // write information of combined SoA buffers back to cell
  combinationSlice[currentSlice]._particleSoABuffer.setViewLength(-1l);
  for (long i = _baseOffsets.front().size() - 2; i >= 0; i--) {
    combinationSlice[currentSlice]._particleSoABuffer.setViewStart(combinationSlicesOffsets[currentSlice][i]);
    combinationSlice[currentSlice]._particleSoABuffer.resizeArrays(combinationSlicesOffsets[currentSlice][i + 1]);
    const long currentOffset = baseIndex + _baseOffsets.front()[i];
    // clear old cell buffer
    cells[currentOffset]._particleSoABuffer.clear();
    // append new cell buffer
    if (combinationSlice[currentSlice]._particleSoABuffer.getNumParticles() != cells[currentOffset].numParticles()) {
      AutoPasLog(error, "1: Count doesn't match. x: {} is: {} should: {}", x,
                 combinationSlice[currentSlice]._particleSoABuffer.getNumParticles(),
                 cells[currentOffset].numParticles());
    }
    cells[currentOffset]._particleSoABuffer.append(combinationSlice[currentSlice]._particleSoABuffer);
  }
  // last cell in stripe need to save whole buffer
  if (x == _cellsPerDimension[0] - _overlap[0] - 1l) {
    for (unsigned long j = 0; j < _overlap[0]; j++) {
      const long sliceNum = (currentSlice + j + 1) % _baseOffsets.size();
      combinationSlice[sliceNum]._particleSoABuffer.setViewLength(-1l);
      for (long i = _baseOffsets[sliceNum].size() - 2; i >= 0; i--) {
        combinationSlice[sliceNum]._particleSoABuffer.setViewStart(combinationSlicesOffsets[sliceNum][i]);
        combinationSlice[sliceNum]._particleSoABuffer.resizeArrays(combinationSlicesOffsets[sliceNum][i + 1]);
        const long currentOffset = baseIndex + _baseOffsets[sliceNum][i];
        // clear old cell buffer
        cells[currentOffset]._particleSoABuffer.clear();
        // append new cell buffer
        if (combinationSlice[sliceNum]._particleSoABuffer.getNumParticles() != cells[currentOffset].numParticles()) {
          AutoPasLog(error, "2: Count doesn't match. x: {} is: {} should: {}", x,
                     combinationSlice[sliceNum]._particleSoABuffer.getNumParticles(),
                     cells[currentOffset].numParticles());
        }
        cells[currentOffset]._particleSoABuffer.append(combinationSlice[sliceNum]._particleSoABuffer);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04CellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::computeOffsets(
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
      for (unsigned long z = 0ul; z <= _overlap[2]; ++z) {
        cellOffsets.push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension));
      }
    }
  }
  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
      _baseOffsets[x].push_back(ov1_squared * x + ov1 * y);
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

  // Create intervals
  const unsigned long numStripes = _cellPairOffsets.size();
  _offsets.resize(numStripes);

  // iterate over all stripes
  for (unsigned long i = 0; i < _cellPairOffsets.size(); ++i) {
    auto &currentStripe = _cellPairOffsets[i];
    // Sort
    std::stable_sort(currentStripe.begin(), currentStripe.end(),
                     [](const auto &a, const auto &b) -> bool { return a.first < b.first; });
    // Collect intervals
    unsigned long current = currentStripe.front().first;
    unsigned long startID = 0;
    for (unsigned long j = 0; j < currentStripe.size(); j++) {
      if (current != currentStripe[j].first) {
        auto interval = std::make_pair(startID, j - 1);
        _offsets[i].push_back(std::make_pair(current, interval));
        startID = j;
      }
    }
    // last interval
    auto interval = std::make_pair(startID, _cellPairOffsets[i].size() - 1);
    _offsets[i].push_back(std::make_pair(current, interval));
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
    _combinationSlicesOffsets.resize(numThreads);
    std::for_each(_combinationSlicesOffsets.begin(), _combinationSlicesOffsets.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _currentSlices.resize(numThreads);
  }
}

}  // namespace autopas
