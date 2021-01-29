/**
 * @file LCC04SoACellHandler.h
 * @author C.Menges
 * @date 04.06.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

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
 * @tparam ParticleCell the type of cells.
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC04SoACellHandler {
 public:
  /**
   * Constructor of the c04 traversal with combined SoA buffers.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   */
  explicit LCC04SoACellHandler(PairwiseFunctor *pairwiseFunctor, std::array<uint64_t, 3> cellsPerDimension,
                               const double interactionLength, const std::array<double, 3> &cellLength,
                               const std::array<uint64_t, 3> &overlap = {1ul, 1ul, 1ul})
      : _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlap(overlap),
        _cellsPerDimension(cellsPerDimension),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(unsigned int)),
        _pairwiseFunctor(pairwiseFunctor) {
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
  void processBaseCell(std::vector<ParticleCell> &cells, uint64_t x, uint64_t y, uint64_t z);

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
  void computeOffsets(std::array<uint64_t, 3> cellsPerDimension);

  /**
   * Interaction Length (cutoff + skin).
   */
  const double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;

  /**
   * overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  const std::array<uint64_t, 3> _overlap;

  const std::array<uint64_t, 3> _cellsPerDimension;

  /**
   * Cells containing combined SoA buffers.
   */
  std::vector<std::vector<ParticleCell>> _combinationSlices;

  /**
   * Current position in circular buffer for each thread.
   */
  std::vector<unsigned int> _currentSlices;

  /**
   * Cell offsets of "base plate".
   */
  std::vector<std::vector<uint64_t>> _baseOffsets;

  /**
   * Type to store an interval of indices.
   */
  using interval_t = std::pair<uint64_t, uint64_t>;

  /**
   * Interactions cell <-> baseplate interval.
   */
  std::vector<std::vector<std::pair<uint64_t, interval_t>>> _offsets;

  /**
   * Partial sums of sizes of combined buffers to determine start and end quickly.
   */
  std::vector<std::vector<std::vector<uint64_t>>> _combinationSlicesOffsets;

  /**
   * Offset factor to avoid false sharing.
   */
  const unsigned int _cacheOffset;

  /**
   * The functor that defines the interaction.
   */
  PairwiseFunctor *_pairwiseFunctor;

  /**
   * Writes buffer content back to cell.
   * @param cells
   * @param baseIndex Index of base cell.
   * @param combinationSlice
   * @param combinationSlicesOffsets
   * @param bufferSlice Index of slice in combinationSlice (source)
   * @param cellSlice Index of slice in _baseOffsets (destination)
   */
  void writeBufferIntoCell(std::vector<ParticleCell> &cells, uint64_t baseIndex,
                           std::vector<ParticleCell> &combinationSlice,
                           std::vector<std::vector<uint64_t>> &combinationSlicesOffsets, uint64_t bufferSlice,
                           uint64_t cellSlice);

  /**
   * Writes cell content into buffer.
   * @param cells
   * @param baseIndex Index of base cell.
   * @param combinationSlice
   * @param combinationSlicesOffsets
   * @param bufferSlice Index of slice in combinationSlice (destination)
   * @param cellSlice Index of slice in _baseOffsets (source)
   */
  void writeCellIntoBuffer(const std::vector<ParticleCell> &cells, uint64_t baseIndex,
                           std::vector<ParticleCell> &combinationSlice,
                           std::vector<std::vector<uint64_t>> &combinationSlicesOffsets, unsigned int bufferSlice,
                           unsigned int cellSlice);

  /**
   * Creates offset intervals (stored in offset) from cellPairOffsets.
   * @param cellPairOffsets Source for interval creation.
   */
  void setupIntervals(std::vector<std::vector<std::pair<uint64_t, uint64_t>>> &cellPairOffsets);
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, uint64_t x, uint64_t y, uint64_t z) {
  const uint64_t baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, _cellsPerDimension);

  // get all information for current thread
  const auto threadID = static_cast<size_t>(autopas_get_thread_num());
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
  for (uint64_t slice = 0; slice < numSlices; slice++) {
    for (auto const &[offset1, interval] : _offsets[(slice + currentSlice) % numSlices]) {
      ParticleCell *cell1 = nullptr;
      size_t cell1ViewStart = 0;
      size_t cell1ViewEnd;

      // special cases (cell1 one is also stored in a combination slice)
      // base cell
      if (offset1 == 0ul) {
        const auto numParticlesBaseCell = cells[baseIndex].numParticles();
        if (numParticlesBaseCell == 0) {
          continue;
        }
        // two cases interval in current  stripe
        cell1 = &combinationSlice[currentSlice];
        auto stripeView = cell1->_particleSoABuffer.constructView(0, numParticlesBaseCell);
        if (slice == currentSlice) {
          // Process stripe with itself
          _pairwiseFunctor->SoAFunctorSingle(stripeView, useNewton3);

          auto restView = cell1->_particleSoABuffer.constructView(numParticlesBaseCell,
                                                                  cell1->_particleSoABuffer.getNumParticles());
          _pairwiseFunctor->SoAFunctorPair(stripeView, restView, useNewton3);
          if (not useNewton3) {
            _pairwiseFunctor->SoAFunctorPair(restView, stripeView, useNewton3);
          }
          cell1ViewEnd = cell1->_particleSoABuffer.getNumParticles();
          continue;
        } else {
          // interval in other stripe
          cell1ViewEnd = numParticlesBaseCell;
        }
      } else if (offset1 == _baseOffsets.front().back()) {
        cell1 = &combinationSlice[currentSlice];
        cell1ViewStart = combinationSlicesOffsets[currentSlice][combinationSlicesOffsets[currentSlice].size() - 2];
        cell1ViewEnd = cell1->_particleSoABuffer.getNumParticles();
      } else if (offset1 == _baseOffsets.back().front()) {
        const auto index = (currentSlice + numSlices - 1) % numSlices;
        if (combinationSlicesOffsets[index][1] == 0) {
          continue;
        }
        cell1 = &combinationSlice[index];
        cell1ViewEnd = combinationSlicesOffsets[index][1];
      } else {
        const uint64_t cellIndex1 = baseIndex + offset1;
        cell1 = &cells[cellIndex1];
        cell1ViewEnd = cell1->_particleSoABuffer.getNumParticles();
      }

      auto &currentCS = combinationSlice[slice];
      const auto &currentCSOffsets = combinationSlicesOffsets[slice];
      auto currentCSViewStart = currentCSOffsets[interval.first];
      auto currentCSViewEnd = currentCSOffsets[interval.second];

      auto cell1View = cell1->_particleSoABuffer.constructView(cell1ViewStart, cell1ViewEnd);
      auto currentCSView = currentCS._particleSoABuffer.constructView(currentCSViewStart, currentCSViewEnd);
      _pairwiseFunctor->SoAFunctorPair(cell1View, currentCSView, useNewton3);
      if (not useNewton3) {
        _pairwiseFunctor->SoAFunctorPair(currentCSView, cell1View, useNewton3);
      }
    }
  }

  // write information of combined SoA buffers back to cell
  writeBufferIntoCell(cells, baseIndex, combinationSlice, combinationSlicesOffsets, currentSlice, 0ul);

  // last cell in stripe need to save whole buffer
  if (x == _cellsPerDimension[0] - _overlap[0] - 1l) {
    for (uint64_t slice = 1; slice < numSlices; slice++) {
      const long bufferSlice = (currentSlice + slice) % numSlices;
      writeBufferIntoCell(cells, baseIndex, combinationSlice, combinationSlicesOffsets, bufferSlice, slice);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::writeCellIntoBuffer(
    const std::vector<ParticleCell> &cells, const uint64_t baseIndex, std::vector<ParticleCell> &combinationSlice,
    std::vector<std::vector<uint64_t>> &combinationSlicesOffsets, const unsigned int bufferSlice,
    const unsigned int cellSlice) {
  // delete old data
  combinationSlice[bufferSlice]._particleSoABuffer.clear();
  combinationSlicesOffsets[bufferSlice].clear();
  // fill buffer and compute partial sums
  uint64_t sum = 0ul;
  combinationSlicesOffsets[bufferSlice].push_back(sum);
  for (const auto offset : _baseOffsets[cellSlice]) {
    const uint64_t otherIndex = baseIndex + offset;
    const ParticleCell &otherCell = cells[otherIndex];
    combinationSlice[bufferSlice]._particleSoABuffer.append(otherCell._particleSoABuffer);
    sum += otherCell.numParticles();
    combinationSlicesOffsets[bufferSlice].push_back(sum);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::writeBufferIntoCell(
    std::vector<ParticleCell> &cells, const uint64_t baseIndex, std::vector<ParticleCell> &combinationSlice,
    std::vector<std::vector<uint64_t>> &combinationSlicesOffsets, const uint64_t bufferSlice,
    const uint64_t cellSlice) {
  auto &buffer = combinationSlice[bufferSlice]._particleSoABuffer;
  for (long i = _baseOffsets[cellSlice].size() - 1l; i >= 0; i--) {
    const auto start = combinationSlicesOffsets[bufferSlice][i];
    const auto end = combinationSlicesOffsets[bufferSlice][i + 1];

    if (start == end) {
      continue;
    }
    buffer.resizeArrays(end);
    auto bufferView = buffer.constructView(start, buffer.getNumParticles());

    const uint64_t currentOffset = baseIndex + _baseOffsets[cellSlice][i];
    // clear old cell buffer
    cells[currentOffset]._particleSoABuffer.clear();
    // make sure everything is correct
    if (bufferView.getNumParticles() != cells[currentOffset].numParticles()) {
      const auto pos = utils::ThreeDimensionalMapping::oneToThreeD(currentOffset, _cellsPerDimension);
      AutoPasLog(error,
                 "Particle number in SoA buffer and cell doesn't match. current position: [{} {} {}] is: {} should: {}",
                 pos[0], pos[1], pos[2], buffer.getNumParticles(), cells[currentOffset].numParticles());
    }
    // append new cell buffer
    cells[currentOffset]._particleSoABuffer.append(bufferView);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::computeOffsets(
    std::array<uint64_t, 3> cellsPerDimension) {
  using std::make_pair;

  std::vector<std::vector<std::pair<uint64_t, uint64_t>>> cellPairOffsets;

  //////////////////////////////
  // @TODO: Replace following lines with vector to support asymmetric cells
  const uint64_t ov1 = _overlap[0] + 1;
  const uint64_t ov1_squared = ov1 * ov1;
  //////////////////////////////

  _baseOffsets.resize(ov1);

  std::array<uint64_t, 3> overlap_1 = utils::ArrayMath::addScalar(_overlap, 1ull);

  std::vector<uint64_t> cellOffsets;
  cellOffsets.reserve(overlap_1[0] * overlap_1[1] * overlap_1[2]);

  cellPairOffsets.clear();
  cellPairOffsets.resize(overlap_1[0]);

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  for (uint64_t x = 0ul; x <= _overlap[0]; ++x) {
    for (uint64_t y = 0ul; y <= _overlap[1]; ++y) {
      _baseOffsets[x].push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, 0ull, cellsPerDimension));
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
          // check whether cell is within cutoff radius
          auto distVec = utils::ArrayMath::mul({std::max(0.0, x - 1.0), std::max(0.0, y - 1.0), std::max(0.0, z - 1.0)},
                                               _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            cellPairOffsets[x].push_back(make_pair(cellOffsets[z], offset));
          }
        }
        // back left
        if (y != _overlap[1] and z != 0) {
          // check whether cell is within cutoff radius
          auto distVec = utils::ArrayMath::mul(
              {std::max(0.0, x - 1.0), std::max(0.0, _overlap[1] - y - 1.0), std::max(0.0, z - 1.0)}, _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            cellPairOffsets[x].push_back(make_pair(cellOffsets[ov1_squared - ov1 + z], offset));
          }
        }
        // front right
        if (x != _overlap[0] and (y != 0 or z != 0)) {
          // check whether cell is within cutoff radius
          auto distVec = utils::ArrayMath::mul(
              {std::max(0.0, _overlap[0] - x - 1.0), std::max(0.0, y - 1.0), std::max(0.0, z - 1.0)}, _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            cellPairOffsets[x].push_back(make_pair(cellOffsets[ov1_squared * _overlap[0] + z], offset));
          }
        }
        // back right
        if (y != _overlap[1] and x != _overlap[0] and z != 0) {
          // check whether cell is within cutoff radius
          auto distVec = utils::ArrayMath::mul(
              {std::max(0.0, _overlap[0] - x - 1.0), std::max(0.0, _overlap[1] - y - 1.0), std::max(0.0, z - 1.0)},
              _cellLength);
          const auto distSquare = utils::ArrayMath::dot(distVec, distVec);
          if (distSquare <= interactionLengthSquare) {
            cellPairOffsets[x].push_back(make_pair(cellOffsets[ov1_squared * ov1 - ov1 + z], offset));
          }
        }
      }
    }
  }
  setupIntervals(cellPairOffsets);
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::setupIntervals(
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> &cellPairOffsets) {
  // Create intervals
  const uint64_t numStripes = cellPairOffsets.size();
  _offsets.resize(numStripes);

  // iterate over all stripes
  for (uint64_t i = 0; i < numStripes; ++i) {
    auto &currentStripe = cellPairOffsets[i];
    // Sort
    std::stable_sort(currentStripe.begin(), currentStripe.end(),
                     [](const auto &a, const auto &b) -> bool { return a.first < b.first; });

    // Collect intervals
    uint64_t current = currentStripe.front().first;
    uint64_t startID =
        std::distance(_baseOffsets[i].begin(),
                      std::find(_baseOffsets[i].begin(), _baseOffsets[i].end(), currentStripe.front().second));
    uint64_t endID = startID;
    for (uint64_t j = 0; j < currentStripe.size(); j++) {
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

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void LCC04SoACellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::resizeBuffers() {
  const auto numThreads = static_cast<size_t>(autopas_get_max_threads());
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
