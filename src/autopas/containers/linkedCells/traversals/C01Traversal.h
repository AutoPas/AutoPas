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
 * newton3 cannot be applied!
 * If combineSoA equals true, SoA buffers are combined slice-wise. Each slice is constructed as FIFO buffer and slices
 * are stored in circular buffer (_combinationSlices).
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 * @tparam combineSoA
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool combineSoA = false>
class C01Traversal
    : public C01BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3, (combineSoA) ? 2 : 3>,
      public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length in CellBlock3D
   */
  explicit C01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                        const double interactionLength = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : C01BasedTraversal < ParticleCell,
      PairwiseFunctor, DataLayout, useNewton3,
      (combineSoA) ? 2 : 3 > (dims, pairwiseFunctor, interactionLength, cellLength),
      _cellFunctor(pairwiseFunctor, interactionLength), _pairwiseFunctor(pairwiseFunctor),
      _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(unsigned int)) {
    computeOffsets();
  }

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  void traverseParticlePairs() override;

  DataLayoutOption getDataLayout() const override { return DataLayout; }

  bool getUseNewton3() const override { return useNewton3; }

  /**
   * C01 traversals are only usable if useNewton3 is disabled and combined SoA buffers are only applicable if SoA is set
   * as DataLayout.
   *
   * This is because the cell functor in the c01 traversal is hardcoded to not allow newton 3 even if only one thread is
   * used.
   *
   * @return
   */
  bool isApplicable() const override {
    return not(DataLayout == DataLayoutOption::cuda) and not useNewton3 and
           not(combineSoA && DataLayout != DataLayoutOption::soa);
  }

  TraversalOption getTraversalType() const override {
    return (combineSoA) ? TraversalOption::c01CombinedSoA : TraversalOption::c01;
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
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, false, false>
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

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3, combineSoA>::computeOffsets() {
  _cellOffsets.resize(2 * this->_overlap[0] + 1);

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  for (long x = -this->_overlap[0]; x <= 0l; ++x) {
    for (long y = -this->_overlap[1]; y <= static_cast<long>(this->_overlap[1]); ++y) {
      for (long z = -this->_overlap[2]; z <= static_cast<long>(this->_overlap[2]); ++z) {
        std::array<double, 3> pos = {};
        pos[0] = std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0];
        pos[1] = std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1];
        pos[2] = std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2];
        const double distSquare = ArrayMath::dot(pos, pos);
        if (distSquare <= interactionLengthSquare) {
          const long currentOffset = utils::ThreeDimensionalMapping::threeToOneD(
              x, y, z, ArrayUtils::static_cast_array<long>(this->_cellsPerDimension));
          const bool containCurrentOffset =
              std::any_of(_cellOffsets[x + this->_overlap[0]].cbegin(), _cellOffsets[x + this->_overlap[0]].cend(),
                          [currentOffset](const auto &e) { return e.first == currentOffset; });
          if (containCurrentOffset) {
            continue;
          }
          for (long ix = x; ix <= std::abs(x); ++ix) {
            const long offset = utils::ThreeDimensionalMapping::threeToOneD(
                ix, y, z, ArrayUtils::static_cast_array<long>(this->_cellsPerDimension));
            const size_t index = ix + this->_overlap[0];
            if (y == 0l and z == 0l) {
              // make sure center of slice is always at the beginning
              _cellOffsets[index].insert(_cellOffsets[index].cbegin(),
                                         std::make_pair(offset, ArrayMath::normalize(pos)));
            } else {
              _cellOffsets[index].push_back(std::make_pair(offset, ArrayMath::normalize(pos)));
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3, combineSoA>::processBaseCell(
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
        _pairwiseFunctor->SoAFunctor(baseCell._particleSoABuffer,
                                     {&(combinationSlice[slice]._particleSoABuffer), startIndex, endIndex}, false);
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

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3, combineSoA>::resizeBuffers() {
  const auto numThreads = static_cast<size_t>(autopas_get_max_threads());
  if (_combinationSlices.size() != numThreads) {
    _combinationSlices.resize(numThreads);
    const auto cellOffsetsSize = _cellOffsets.size();
    std::for_each(_combinationSlices.begin(), _combinationSlices.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _currentSlices.resize(numThreads * _cacheOffset);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool combineSoA>
inline void C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3, combineSoA>::traverseParticlePairs() {
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
