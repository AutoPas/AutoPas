/**
 * @file C01Traversal.h
 * @author nguyen
 * @date 16.09.2018
 */

#pragma once

#include "LinkedCellTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c01 traversal and the c01 traversal with combined SoA buffers.
 *
 * The traversal uses the c01 base step performed on every single cell.
 * \image html C01.png "C01 base step in 2D. (dark blue cell = base cell)"
 * newton3 cannot be applied!
 * If combineSoA equals true, SoA buffers are combined slice-wise, as described below.
 *
 * Each slice is constructed as FIFO buffer and slices
 * are stored in circular buffer (_combinationSlices).
 *
 * <b>Combination Principle</b>
 *
 * The sphere of interacting cells is divided into slices along the y-axis, as seen in Figure 1. Each
slice represents a combined SoA buffer. We assume that the creation of combined SoA buffers starts at the
first evaluated cell of an axis. Since there is no previously evaluated cell on the beginning of the axis, the whole
buffer is initially filled by copying the data from the cells. The combined SoA buffers are stored in a circular/ring
buffer (_combinationSlices), shown in Figure 2. To keep track of the position of the first slice inside the circular
buffer, an additional variable for the start index (_currentSlices) is necessary which is initialized to zero. Inside
the combined SoA buffer, the particles are ordered according their insertion. Note, that the base cell (dark blue)
always represents the first cell in the slice buffer.

Now, all interactions can be calculated by iterating over all slices and compute the interactions with the current base
cell. Since information in the buffers is not persistent, it is important to use the SoA buffer of the base cell
and not the copy inside of the combined buffer. If the current base cell interacts with the slice which contains a copy
of the base cell, it is necessary to exclude the copy from the calculations. Otherwise, particles would interact with
themselves. It is not possible to remove the copy from the buffer slice, since the current base cell is an interacting
cell of the next base cell. Therefore, the SoA data structure defines a custom view on
the underlying data structure. Since the copy of the current base cell is the first cell in the buffer slice, it is
possible to set the start of the view to the first particle after the base cell.

In the next step, the sphere of interactions moves one cell further. Figure 1 shows that most of the new
and old interaction cells are the same. To reduce the number of copies, we'll keep the combined SoA buffers from
the previous step and only apply an update to them. Here, the symmetry of interactions can be exploited. Each combined
SoA buffer is a LIFO (Last In, First Out) buffer, meaning that they are demolished in the reverse order of their
construction. This effect can be seen as well in the circular buffer (Figure 2). The former leftmost slice goes
completely out of scope and can be deleted from the circular buffer. The position of this slice inside the circular
buffer is represented by the start index inside the buffer. At the same time, a new SoA buffer is created to represent
the rightmost slice, which has just entered the interaction sphere. This slice is written on the same index inside the
SoA buffer as the former leftmost slice. Since the first slice has moved inside the circular buffer, the start index is
incremented. Due to the circular characteristics, an additional modulo operation is applied to the start index to jump
back to the first slice in the ring buffer if the end is reached. At this point, the buffer is fully updated and all
interactions can be evaluated. This procedure is repeated until the end of the axis is reached.

Since the offsets of the interacting cells relative to the base cell do not change during the traversal, they can be
computed beforehand. All offsets are stored in a 2D-array (_cellOffsets) where the first dimension represents the
individual slices of the interaction sphere and the second dimension represents the cell offsets inside the slice. The
cell offsets are sorted to resemble the order of growth/destruction of the combined SoA buffers. This is important since
the initialized buffer must show the same behavior as a buffer which was updated multiple times.
 *
 * \image html C01_combined.png "Figure 1: Movement of the interaction sphere in 2D. (dark blue cell = base cell)"
 * \image html C01_combined_cache.png "Figure 2: Evolution of the circular buffer."
 * \image html C01_combined_legend.png
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 * @tparam combineSoA
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          bool combineSoA = false>
class C01Traversal
    : public C01BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, (combineSoA ? 2 : 3)>,
      public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length in CellBlock3D
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit C01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                        const double interactionLength, const std::array<double, 3> &cellLength)
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, (combineSoA ? 2 : 3)>(
            dims, pairwiseFunctor, interactionLength, cellLength),
        _cellFunctor(pairwiseFunctor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/),
        _pairwiseFunctor(pairwiseFunctor),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(unsigned int)) {
    computeOffsets();
  }

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  void traverseParticlePairs() override;

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  /**
   * C01 traversals are only usable if useNewton3 is disabled and combined SoA buffers are only applicable if SoA is set
   * as DataLayout.
   *
   * This is because the cell functor in the c01 traversal is hardcoded to not allow newton 3 even if only one thread is
   * used.
   *
   * @return
   */
  [[nodiscard]] bool isApplicable() const override {
    return not(dataLayout == DataLayoutOption::cuda) and not useNewton3 and
           not(combineSoA && dataLayout != DataLayoutOption::soa);
  }

  [[nodiscard]] TraversalOption getTraversalType() const override {
    return (combineSoA) ? TraversalOption::lc_c01_combined_SoA : TraversalOption::lc_c01;
  }

 private:
  /**
   * Computes all interactions between the base
   * cell and adjacent cells.
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  inline void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  /**
   * Appends all needed Attributes to the SoA buffer in cell.
   * @tparam I
   * @param cell
   * @param appendCell
   */
  template <std::size_t... I>
  inline constexpr void appendNeeded(ParticleCell &cell, ParticleCell &appendCell, std::index_sequence<I...>) {
    cell._particleSoABuffer.template append<std::get<I>(PairwiseFunctor::getNeededAttr(std::false_type()))...>(
        appendCell._particleSoABuffer);
  }

  /**
   * Resizes all buffers needed for combined SoA buffers (_combinationSlices, _currentSlices) to fit the current number
   * of threads and cell offsets (= _cellOffsets.size()).
   */
  void resizeBuffers();

  /**
   * Pairs for processBaseCell().
   * @note std::map not applicable since ordering arising from insertion is important for later processing!
   */
  std::vector<std::vector<std::pair<long, std::array<double, 3>>>> _cellOffsets;

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, dataLayout, false, false>
      _cellFunctor;

  PairwiseFunctor *_pairwiseFunctor;

  /**
   * Cells containing combined SoA buffers.
   */
  std::vector<std::vector<ParticleCell>> _combinationSlices;

  /**
   * Current index in _combinationSlices.
   */
  std::vector<unsigned int> _currentSlices;

  /**
   * Offset factor to avoid false sharing.
   */
  const unsigned int _cacheOffset;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, combineSoA>::computeOffsets() {
  _cellOffsets.resize(2 * this->_overlap[0] + 1);

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  for (long x = -this->_overlap[0]; x <= 0l; ++x) {
    for (long y = -this->_overlap[1]; y <= static_cast<long>(this->_overlap[1]); ++y) {
      for (long z = -this->_overlap[2]; z <= static_cast<long>(this->_overlap[2]); ++z) {
        std::array<double, 3> pos = {};
        pos[0] = std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0];
        pos[1] = std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1];
        pos[2] = std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2];
        const double distSquare = utils::ArrayMath::dot(pos, pos);
        if (distSquare <= interactionLengthSquare) {
          const long currentOffset = utils::ThreeDimensionalMapping::threeToOneD(
              x, y, z, utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension));
          const bool containCurrentOffset =
              std::any_of(_cellOffsets[x + this->_overlap[0]].cbegin(), _cellOffsets[x + this->_overlap[0]].cend(),
                          [currentOffset](const auto &e) { return e.first == currentOffset; });
          if (containCurrentOffset) {
            continue;
          }
          for (long ix = x; ix <= std::abs(x); ++ix) {
            const long offset = utils::ThreeDimensionalMapping::threeToOneD(
                ix, y, z, utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension));
            const size_t index = ix + this->_overlap[0];
            if (y == 0l and z == 0l) {
              // make sure center of slice is always at the beginning
              _cellOffsets[index].insert(_cellOffsets[index].cbegin(),
                                         std::make_pair(offset, utils::ArrayMath::normalize(pos)));
            } else {
              _cellOffsets[index].push_back(std::make_pair(offset, utils::ArrayMath::normalize(pos)));
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, combineSoA>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
  ParticleCell &baseCell = cells[baseIndex];
  const size_t cOffSize = _cellOffsets.size();

  if constexpr (combineSoA) {
    // Iteration along x

    const auto threadID = static_cast<size_t>(autopas_get_thread_num());
    auto &currentSlice = _currentSlices[threadID * _cacheOffset];
    auto &combinationSlice = _combinationSlices[threadID];

    // First cell needs to initialize whole buffer
    if (x == this->_overlap[0]) {
      currentSlice = 0;
      for (unsigned int offsetSlice = 0; offsetSlice < cOffSize; offsetSlice++) {
        combinationSlice[offsetSlice]._particleSoABuffer.clear();
        for (const auto &offset : _cellOffsets[offsetSlice]) {
          const unsigned long otherIndex = baseIndex + offset.first;
          ParticleCell &otherCell = cells[otherIndex];
          appendNeeded(combinationSlice[offsetSlice], otherCell,
                       std::make_index_sequence<PairwiseFunctor::getNeededAttr(std::false_type()).size()>{});
        }
      }
    } else {
      // reduce size
      size_t i = 0;
      const size_t midSlice = (currentSlice + this->_overlap[0] + 1) % cOffSize;
      for (size_t slice = (currentSlice + 1) % cOffSize; slice != midSlice; ++slice %= cOffSize, ++i) {
        size_t newSize = 0;
        for (const auto &offset : _cellOffsets[i]) {
          const unsigned long otherIndex = baseIndex + offset.first;
          ParticleCell &otherCell = cells[otherIndex];
          newSize += otherCell.numParticles();
        }
        combinationSlice[slice]._particleSoABuffer.resizeArrays(newSize);
      }
      // append buffers
      for (size_t slice = midSlice; slice != currentSlice; ++slice %= cOffSize, ++i) {
        for (auto offsetIndex = _cellOffsets[(i + 1) % cOffSize].size(); offsetIndex < _cellOffsets[i].size();
             ++offsetIndex) {
          const unsigned long otherIndex = baseIndex + _cellOffsets[i][offsetIndex].first;
          ParticleCell &otherCell = cells[otherIndex];
          appendNeeded(combinationSlice[slice], otherCell,
                       std::make_index_sequence<PairwiseFunctor::getNeededAttr(std::false_type()).size()>{});
        }
      }

      combinationSlice[currentSlice]._particleSoABuffer.clear();

      for (const auto &offset : _cellOffsets.back()) {
        const unsigned long otherIndex = baseIndex + offset.first;
        ParticleCell &otherCell = cells[otherIndex];
        appendNeeded(combinationSlice[currentSlice], otherCell,
                     std::make_index_sequence<PairwiseFunctor::getNeededAttr(std::false_type()).size()>{});
      }

      ++currentSlice %= cOffSize;
    }

    // calculate all interactions
    for (unsigned int slice = 0; slice < cOffSize; slice++) {
      if (slice == (currentSlice + this->_overlap[0]) % cOffSize) {
        // slice contains base cell -> skip particles of base cell. This is not supported by CellFunctor, so call
        // pairwise functor directly.
        auto startIndex = baseCell.numParticles();
        auto endIndex = combinationSlice[slice]._particleSoABuffer.getNumParticles();
        _pairwiseFunctor->SoAFunctorPair(baseCell._particleSoABuffer,
                                         {&(combinationSlice[slice]._particleSoABuffer), startIndex, endIndex}, false,
                                         true);
        // compute base cell
        this->_cellFunctor.processCell(baseCell);
      } else {
        this->_cellFunctor.processCellPair(baseCell, combinationSlice[slice]);
      }
    }
  } else {
    for (const auto &slice : _cellOffsets) {
      for (auto const &[offset, r] : slice) {
        const unsigned long otherIndex = baseIndex + offset;
        ParticleCell &otherCell = cells[otherIndex];

        if (baseIndex == otherIndex) {
          this->_cellFunctor.processCell(baseCell);
        } else {
          this->_cellFunctor.processCellPair(baseCell, otherCell, r);
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, combineSoA>::resizeBuffers() {
  const auto numThreads = static_cast<size_t>(autopas_get_max_threads());
  if (_combinationSlices.size() != numThreads) {
    _combinationSlices.resize(numThreads);
    const auto cellOffsetsSize = _cellOffsets.size();
    std::for_each(_combinationSlices.begin(), _combinationSlices.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _currentSlices.resize(numThreads * _cacheOffset);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, combineSoA>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  if (not this->isApplicable()) {
    if constexpr (combineSoA) {
      utils::ExceptionHandler::exception(
          "The C01 traversal with combined SoA buffers cannot work with data layout AoS and enabled newton3 (unless "
          "only one thread is used)!");
    } else {
      utils::ExceptionHandler::exception(
          "The C01 traversal cannot work with enabled newton3 (unless only one thread is used)!");
    }
  }
  if constexpr (combineSoA) {
    resizeBuffers();
  }
  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
