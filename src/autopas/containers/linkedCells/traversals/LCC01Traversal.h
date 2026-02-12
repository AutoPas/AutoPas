/**
 * @file LCC01Traversal.h
 * @author nguyen
 * @date 16.09.2018
 */

#pragma once

#include "LCTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/containers/cellTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/checkFunctorType.h"

namespace autopas {

/**
 * This class provides the c01 traversal and the c01 traversal with combined SoA buffers.
 *
 * The traversal uses the c01 base step performed on every single cell.
 * \image html C01.png "C01 base step in 2D. (dark blue cell = base cell)"
 * newton3 cannot be applied!
 * If combineSoA equals true, SoA buffers are combined slice-wise, as described below.
 * Note: combineSoA is only implemented for pairwise Functors at the moment!
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
 * @tparam Functor The functor that defines the interaction of two particles.
 * @tparam combineSoA
 */
template <class ParticleCell, class Functor, bool combineSoA = false>
class LCC01Traversal : public C01BasedTraversal<ParticleCell, Functor, (combineSoA ? 2 : 3)>,
                       public LCTraversalInterface {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param functor The functor that defines the interaction of particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length in CellBlock3D
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC01Traversal(const std::array<unsigned long, 3> &dims, Functor *functor, const double interactionLength,
                          const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : C01BasedTraversal<ParticleCell, Functor, (combineSoA ? 2 : 3)>(dims, functor, interactionLength, cellLength,
                                                                       dataLayout, useNewton3),
        _cellFunctor(functor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/,
                     dataLayout, useNewton3),
        _functor(functor),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(unsigned int)) {
    this->computeOffsets();
  }

  void traverseParticles() override;

  /**
   * C01 traversals are only usable if useNewton3 is disabled and combined SoA buffers are only applicable if SoA is set
   * as DataLayout.
   *
   * This is because the cell functor in the c01 traversal is hardcoded to not allow newton 3 even if only one thread is
   * used.
   *
   * Also, combined SoA buffers are only implemented for pairwise interactions.
   *
   * @return
   */
  [[nodiscard]] bool isApplicable() const override {
    return not this->_useNewton3 and not(combineSoA and this->_dataLayout != DataLayoutOption::soa) and
           not(combineSoA and not utils::isPairwiseFunctor<Functor>());
  }

  [[nodiscard]] TraversalOption getTraversalType() const override {
    return (combineSoA) ? TraversalOption::lc_c01_combined_SoA : TraversalOption::lc_c01;
  }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellFunctor.setSortingThreshold(sortingThreshold); }

 private:
  // CellFunctor type for either Pairwise or Triwise Functors.
  using CellFunctorType = std::conditional_t<decltype(utils::isPairwiseFunctor<Functor>())::value,
                                             internal::CellFunctor<ParticleCell, Functor, /*bidirectional*/ false>,
                                             internal::CellFunctor3B<ParticleCell, Functor, /*bidirectional*/ false>>;

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
   * Pairwise implementation of processBaseCell().
   * @copydoc processBaseCell()
   */
  inline void processBaseCellPairwise(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y,
                                      unsigned long z);

  /**
   * Triwise implementation of processBaseCell().
   * @copydoc processBaseCell()
   */
  inline void processBaseCellTriwise(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y,
                                     unsigned long z);

  /**
   * Appends all needed Attributes to the SoA buffer in cell.
   * @tparam I
   * @param cell
   * @param appendCell
   */
  template <std::size_t... I>
  constexpr void appendNeeded(ParticleCell &cell, ParticleCell &appendCell, std::index_sequence<I...>) {
    cell._particleSoABuffer.template append<std::get<I>(Functor::getNeededAttr(std::false_type()))...>(
        appendCell._particleSoABuffer);
  }

  /**
   * Resizes all buffers needed for combined SoA buffers (_combinationSlices, _currentSlices) to fit the current number
   * of threads and cell offsets (= _cellOffsets.size()).
   */
  void resizeBuffers();

  /**
   * CellFunctor to be used for the traversal defining the interaction between two or more cells.
   */
  CellFunctorType _cellFunctor;

  /**
   * Functor defining pairwise or triwise particle interactions.
   */
  Functor *_functor;

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

template <class ParticleCell, class Functor, bool combineSoA>
inline void LCC01Traversal<ParticleCell, Functor, combineSoA>::processBaseCell(std::vector<ParticleCell> &cells,
                                                                               unsigned long x, unsigned long y,
                                                                               unsigned long z) {
  if constexpr (utils::isPairwiseFunctor<Functor>()) {
    processBaseCellPairwise(cells, x, y, z);
  } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
    processBaseCellTriwise(cells, x, y, z);
  } else {
    utils::ExceptionHandler::exception(
        "LCC01Traversal::processBaseCell(): Functor {} is not of type PairwiseFunctor or TriwiseFunctor.",
        _functor->getName());
  }
}

template <class ParticleCell, class Functor, bool combineSoA>
inline void LCC01Traversal<ParticleCell, Functor, combineSoA>::processBaseCellPairwise(std::vector<ParticleCell> &cells,
                                                                                       unsigned long x, unsigned long y,
                                                                                       unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
  ParticleCell &baseCell = cells[baseIndex];
  const size_t cOffSize = this->_cellOffsets.size();

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
        for (const auto &offset : this->_cellOffsets[offsetSlice]) {
          const unsigned long otherIndex = baseIndex + offset.first;
          ParticleCell &otherCell = cells[otherIndex];
          appendNeeded(combinationSlice[offsetSlice], otherCell,
                       std::make_index_sequence<Functor::getNeededAttr(std::false_type()).size()>{});
        }
      }
    } else {
      // reduce size
      size_t i = 0;
      const size_t midSlice = (currentSlice + this->_overlap[0] + 1) % cOffSize;
      for (size_t slice = (currentSlice + 1) % cOffSize; slice != midSlice; ++slice %= cOffSize, ++i) {
        size_t newSize = 0;
        for (const auto &offset : this->_cellOffsets[i]) {
          const unsigned long otherIndex = baseIndex + offset.first;
          ParticleCell &otherCell = cells[otherIndex];
          newSize += otherCell.size();
        }
        combinationSlice[slice]._particleSoABuffer.resizeArrays(newSize);
      }
      // append buffers
      for (size_t slice = midSlice; slice != currentSlice; ++slice %= cOffSize, ++i) {
        for (auto offsetIndex = this->_cellOffsets[(i + 1) % cOffSize].size(); offsetIndex < this->_cellOffsets[i].size();
             ++offsetIndex) {
          const unsigned long otherIndex = baseIndex + this->_cellOffsets[i][offsetIndex].first;
          ParticleCell &otherCell = cells[otherIndex];
          appendNeeded(combinationSlice[slice], otherCell,
                       std::make_index_sequence<Functor::getNeededAttr(std::false_type()).size()>{});
        }
      }

      combinationSlice[currentSlice]._particleSoABuffer.clear();

      for (const auto &offset : this->_cellOffsets.back()) {
        const unsigned long otherIndex = baseIndex + offset.first;
        ParticleCell &otherCell = cells[otherIndex];
        appendNeeded(combinationSlice[currentSlice], otherCell,
                     std::make_index_sequence<Functor::getNeededAttr(std::false_type()).size()>{});
      }

      ++currentSlice %= cOffSize;
    }

    // calculate all interactions
    for (unsigned int slice = 0; slice < cOffSize; slice++) {
      if (slice == (currentSlice + this->_overlap[0]) % cOffSize) {
        // slice contains base cell -> skip particles of base cell. This is not supported by CellFunctor, so call
        // pairwise functor directly.
        auto startIndex = baseCell.size();
        auto endIndex = combinationSlice[slice]._particleSoABuffer.size();
        _functor->SoAFunctorPair(baseCell._particleSoABuffer,
                                 {&(combinationSlice[slice]._particleSoABuffer), startIndex, endIndex}, false);
        // compute base cell
        this->_cellFunctor.processCell(baseCell);
      } else {
        this->_cellFunctor.processCellPair(baseCell, combinationSlice[slice]);
      }
    }
  } else {
    for (const auto &slice : this->_cellOffsets) {
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

template <class ParticleCell, class Functor, bool combineSoA>
inline void LCC01Traversal<ParticleCell, Functor, combineSoA>::processBaseCellTriwise(std::vector<ParticleCell> &cells,
                                                                                      unsigned long x, unsigned long y,
                                                                                      unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
  ParticleCell &baseCell = cells[baseIndex];

  for (auto const &[offset1, offset2, r] : this->_cellOffsets) {
    const unsigned long otherIndex1 = baseIndex + offset1;
    const unsigned long otherIndex2 = baseIndex + offset2;
    ParticleCell &otherCell1 = cells[otherIndex1];
    ParticleCell &otherCell2 = cells[otherIndex2];

    if (baseIndex == otherIndex1 and baseIndex == otherIndex2) {
      this->_cellFunctor.processCell(baseCell);
    } else if (baseIndex == otherIndex1 and baseIndex != otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherCell2);
    } else if (baseIndex != otherIndex1 and baseIndex == otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherCell1);
    } else if (baseIndex != otherIndex1 and otherIndex1 == otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherCell1);
    } else {
      this->_cellFunctor.processCellTriple(baseCell, otherCell1, otherCell2, r);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool combineSoA>
inline void LCC01Traversal<ParticleCell, PairwiseFunctor, combineSoA>::resizeBuffers() {
  const auto numThreads = static_cast<size_t>(autopas_get_max_threads());
  if (_combinationSlices.size() != numThreads) {
    _combinationSlices.resize(numThreads);
    const auto cellOffsetsSize = this->_cellOffsets.size();
    std::for_each(_combinationSlices.begin(), _combinationSlices.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _currentSlices.resize(numThreads * _cacheOffset);
  }
}

template <class ParticleCell, class Functor, bool combineSoA>
inline void LCC01Traversal<ParticleCell, Functor, combineSoA>::traverseParticles() {
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
