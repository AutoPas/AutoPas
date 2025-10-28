/**
 * @file LCC04SoACellHandler.h
 * @author C.Menges
 * @date 04.06.2019
 */

#pragma once

#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandlerUtility.h"
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
 */
template <class ParticleCell, class PairwiseFunctor>
class LCC04SoACellHandler {
 public:
  /**
   * Constructor of the c04 traversal with combined SoA buffers.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit LCC04SoACellHandler(PairwiseFunctor *pairwiseFunctor, const std::array<unsigned long, 3> &cellsPerDimension,
                               double interactionLength, const std::array<double, 3> &cellLength,
                               DataLayoutOption dataLayout, bool useNewton3,
                               const std::array<unsigned long, 3> &overlap = {1ul, 1ul, 1ul})
      : _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlap(overlap),
        _cellsPerDimension(cellsPerDimension),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(unsigned int)),
        _pairwiseFunctor(pairwiseFunctor),
        _dataLayout(dataLayout),
        _useNewton3(useNewton3) {
    setupIntervals(cellsPerDimension);
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
  const std::array<unsigned long, 3> _overlap;

  const std::array<unsigned long, 3> _cellsPerDimension;

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
  std::vector<std::vector<unsigned long>> _baseOffsets;

  /**
   * Type to store an interval of indices.
   */
  using interval_t = std::pair<unsigned long, unsigned long>;

  /**
   * Interactions cell <-> baseplate interval.
   */
  std::vector<std::vector<std::pair<unsigned long, interval_t>>> _offsets;

  /**
   * Partial sums of sizes of combined buffers to determine start and end quickly.
   */
  std::vector<std::vector<std::vector<unsigned long>>> _combinationSlicesOffsets;

  /**
   * Offset factor to avoid false sharing.
   */
  const unsigned int _cacheOffset;

  /**
   * The functor that defines the interaction.
   */
  PairwiseFunctor *_pairwiseFunctor;

  /**
   * The datalayout used by this traversal.
   */
  DataLayoutOption _dataLayout;

  /**
   * If this traversal makes use of newton3.
   */
  bool _useNewton3;

  /**
   * Writes buffer content back to cell.
   * @param cells
   * @param baseIndex Index of base cell.
   * @param combinationSlice
   * @param combinationSlicesOffsets
   * @param bufferSlice Index of slice in combinationSlice (source)
   * @param cellSlice Index of slice in _baseOffsets (destination)
   */
  void writeBufferIntoCell(std::vector<ParticleCell> &cells, unsigned long baseIndex,
                           std::vector<ParticleCell> &combinationSlice,
                           std::vector<std::vector<unsigned long>> &combinationSlicesOffsets, unsigned long bufferSlice,
                           unsigned long cellSlice);

  /**
   * Writes cell content into buffer.
   * @param cells
   * @param baseIndex Index of base cell.
   * @param combinationSlice
   * @param combinationSlicesOffsets
   * @param bufferSlice Index of slice in combinationSlice (destination)
   * @param cellSlice Index of slice in _baseOffsets (source)
   */
  void writeCellIntoBuffer(const std::vector<ParticleCell> &cells, unsigned long baseIndex,
                           std::vector<ParticleCell> &combinationSlice,
                           std::vector<std::vector<unsigned long>> &combinationSlicesOffsets, unsigned int bufferSlice,
                           unsigned int cellSlice);

  /**
   * Creates offset intervals (stored in offset) from cellPairOffsets.
   * @param cellsPerDimension cells per dimension.
   */
  void setupIntervals(const std::array<unsigned long, 3> &cellsPerDimension);
};

template <class ParticleCell, class PairwiseFunctor>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor>::processBaseCell(std::vector<ParticleCell> &cells,
                                                                                unsigned long x, unsigned long y,
                                                                                unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, _cellsPerDimension);

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
  for (unsigned long slice = 0; slice < numSlices; slice++) {
    for (auto const &[offset1, interval] : _offsets[(slice + currentSlice) % numSlices]) {
      ParticleCell *cell1 = nullptr;
      size_t cell1ViewStart = 0;
      size_t cell1ViewEnd;

      // special cases (cell1 one is also stored in a combination slice)
      // base cell
      if (offset1 == 0ul) {
        const auto numParticlesBaseCell = cells[baseIndex].size();
        if (numParticlesBaseCell == 0) {
          continue;
        }
        // two cases interval in current  stripe
        cell1 = &combinationSlice[currentSlice];
        auto stripeView = cell1->_particleSoABuffer.constructView(0, numParticlesBaseCell);
        if (slice == currentSlice) {
          // Process stripe with itself
          _pairwiseFunctor->SoAFunctorSingle(stripeView, _useNewton3);

          auto restView =
              cell1->_particleSoABuffer.constructView(numParticlesBaseCell, cell1->_particleSoABuffer.size());
          _pairwiseFunctor->SoAFunctorPair(stripeView, restView, _useNewton3);
          if (not _useNewton3) {
            _pairwiseFunctor->SoAFunctorPair(restView, stripeView, _useNewton3);
          }
          cell1ViewEnd = cell1->_particleSoABuffer.size();
          continue;
        } else {
          // interval in other stripe
          cell1ViewEnd = numParticlesBaseCell;
        }
      } else if (offset1 == _baseOffsets.front().back()) {
        cell1 = &combinationSlice[currentSlice];
        cell1ViewStart = combinationSlicesOffsets[currentSlice][combinationSlicesOffsets[currentSlice].size() - 2];
        cell1ViewEnd = cell1->_particleSoABuffer.size();
      } else if (offset1 == _baseOffsets.back().front()) {
        const auto index = (currentSlice + numSlices - 1) % numSlices;
        if (combinationSlicesOffsets[index][1] == 0) {
          continue;
        }
        cell1 = &combinationSlice[index];
        cell1ViewEnd = combinationSlicesOffsets[index][1];
      } else {
        const unsigned long cellIndex1 = baseIndex + offset1;
        cell1 = &cells[cellIndex1];
        cell1ViewEnd = cell1->_particleSoABuffer.size();
      }

      auto &currentCS = combinationSlice[slice];
      const auto &currentCSOffsets = combinationSlicesOffsets[slice];
      auto currentCSViewStart = currentCSOffsets[interval.first];
      auto currentCSViewEnd = currentCSOffsets[interval.second];

      auto cell1View = cell1->_particleSoABuffer.constructView(cell1ViewStart, cell1ViewEnd);
      auto currentCSView = currentCS._particleSoABuffer.constructView(currentCSViewStart, currentCSViewEnd);
      _pairwiseFunctor->SoAFunctorPair(cell1View, currentCSView, _useNewton3);
      if (not _useNewton3) {
        _pairwiseFunctor->SoAFunctorPair(currentCSView, cell1View, _useNewton3);
      }
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

template <class ParticleCell, class PairwiseFunctor>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor>::writeCellIntoBuffer(
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
    sum += otherCell.size();
    combinationSlicesOffsets[bufferSlice].push_back(sum);
  }
}

template <class ParticleCell, class PairwiseFunctor>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor>::writeBufferIntoCell(
    std::vector<ParticleCell> &cells, const unsigned long baseIndex, std::vector<ParticleCell> &combinationSlice,
    std::vector<std::vector<unsigned long>> &combinationSlicesOffsets, const unsigned long bufferSlice,
    const unsigned long cellSlice) {
  auto &buffer = combinationSlice[bufferSlice]._particleSoABuffer;
  for (long i = _baseOffsets[cellSlice].size() - 1l; i >= 0; i--) {
    const auto start = combinationSlicesOffsets[bufferSlice][i];
    const auto end = combinationSlicesOffsets[bufferSlice][i + 1];

    if (start == end) {
      continue;
    }
    buffer.resizeArrays(end);
    auto bufferView = buffer.constructView(start, buffer.size());

    const unsigned long currentOffset = baseIndex + _baseOffsets[cellSlice][i];
    // clear old cell buffer
    cells[currentOffset]._particleSoABuffer.clear();
    // make sure everything is correct
    if (bufferView.size() != cells[currentOffset].size()) {
      const auto pos = utils::ThreeDimensionalMapping::oneToThreeD(currentOffset, _cellsPerDimension);
      AutoPasLog(ERROR,
                 "Particle number in SoA buffer and cell doesn't match. current position: [{} {} {}] is: {} should: {}",
                 pos[0], pos[1], pos[2], buffer.size(), cells[currentOffset].size());
    }
    // append new cell buffer
    cells[currentOffset]._particleSoABuffer.append(bufferView);
  }
}

template <class ParticleCell, class PairwiseFunctor>
inline void LCC04SoACellHandler<ParticleCell, PairwiseFunctor>::setupIntervals(
    const std::array<unsigned long, 3> &cellsPerDimension) {
  _baseOffsets.resize(_overlap[0] + 1);
  for (unsigned long x = 0ul; x <= _overlap[0]; ++x) {
    for (unsigned long y = 0ul; y <= _overlap[1]; ++y) {
      _baseOffsets[x].push_back(utils::ThreeDimensionalMapping::threeToOneD(x, y, 0ul, cellsPerDimension));
    }
  }

  std::vector<LCC08CellHandlerUtility::OffsetPairVector> cellPairOffsets =
      LCC08CellHandlerUtility::computePairwiseCellOffsetsC08<LCC08CellHandlerUtility::C08OffsetMode::c04CellPairs>(
          cellsPerDimension, this->_cellLength, this->_interactionLength);

  // Create intervals
  const unsigned long numStripes = cellPairOffsets.size();
  _offsets.resize(numStripes);

  // iterate over all stripes
  for (unsigned long i = 0; i < numStripes; ++i) {
    auto &currentStripe = cellPairOffsets[i];
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

template <class ParticleCell, class PairwiseFunctor>
void LCC04SoACellHandler<ParticleCell, PairwiseFunctor>::resizeBuffers() {
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
