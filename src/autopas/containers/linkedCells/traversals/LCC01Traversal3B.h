/**
 * @file LCC01Traversal3B.h
 * @author muehlhaeusser
 * @date 05.09.2023
 */

#pragma once

#include "autopas/containers/TriwiseTraversalInterface.h"
#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor3B.h"
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
 * @tparam Functor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC01Traversal3B
    : public C01BasedTraversal<ParticleCell, Functor, dataLayout, useNewton3, 3>,
      public LCTraversalInterface<ParticleCell>, public TriwiseTraversalInterface {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param functor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length in CellBlock3D
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC01Traversal3B(const std::array<unsigned long, 3> &dims, Functor *functor,
                          const double interactionLength, const std::array<double, 3> &cellLength)
      : C01BasedTraversal<ParticleCell, Functor, dataLayout, useNewton3, 3>(
            dims, functor, interactionLength, cellLength),
        _cellFunctor(functor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/),
        _functor(functor),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(unsigned int)) {
    computeOffsets();
  }

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  void traverseParticleTriplets() override;

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
    return not useNewton3;
  }

  [[nodiscard]] TraversalOption getTraversalType() const override {
    return TraversalOption::lc_c01_3b;
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
    cell._particleSoABuffer.template append<std::get<I>(Functor::getNeededAttr(std::false_type()))...>(
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
  std::vector<std::tuple<long, long, std::array<double, 3>>> _cellOffsets;

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor3B<typename ParticleCell::ParticleType, ParticleCell, Functor, dataLayout, false, false>
      _cellFunctor;

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

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC01Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::computeOffsets() {
  using namespace utils::ArrayMath::literals;

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

  // offsets for the first cell
  for (long x1 = -this->_overlap[0]; x1 <= static_cast<long>(this->_overlap[0]); ++x1) {
    for (long y1 = -this->_overlap[1]; y1 <= static_cast<long>(this->_overlap[1]); ++y1) {
      for (long z1 = -this->_overlap[2]; z1 <= static_cast<long>(this->_overlap[2]); ++z1) {
        // check distance between base cell and cell 1
        const std::array<double, 3> dist01 = {
            std::max(0l, (std::abs(x1) - 1l)) * this->_cellLength[0],
            std::max(0l, (std::abs(y1) - 1l)) * this->_cellLength[1],
            std::max(0l, (std::abs(z1) - 1l)) * this->_cellLength[2],
        };
        const double distSquare = utils::ArrayMath::dot(dist01, dist01);
        if (distSquare <= interactionLengthSquare) break;

        // offsets for the second cell
        for (long x2 = x1; x2 <= static_cast<long>(this->_overlap[0]); ++x2) {
          for (long y2 = y1; y2 <= static_cast<long>(this->_overlap[1]); ++y2) {
            for (long z2 = z1; z2 <= static_cast<long>(this->_overlap[2]); ++z2) {
              // check distance between cell 1 and cell 2
              const std::array<double, 3> dist12 = {
                  std::max(0l, (std::abs(x1 - x2) - 1l)) * this->_cellLength[0],
                  std::max(0l, (std::abs(y1 - y2) - 1l)) * this->_cellLength[1],
                  std::max(0l, (std::abs(z1 - z2) - 1l)) * this->_cellLength[2],
              };
              const double dist12Squared = utils::ArrayMath::dot(dist12, dist12);
              if (dist12Squared <= interactionLengthSquare) break;

              // check distance between base cell and cell 2
              const std::array<double, 3> dist02 = {
                  std::max(0l, (std::abs(x2) - 1l)) * this->_cellLength[0],
                  std::max(0l, (std::abs(y2) - 1l)) * this->_cellLength[1],
                  std::max(0l, (std::abs(z2) - 1l)) * this->_cellLength[2],
              };
              const double dist02Squared = utils::ArrayMath::dot(dist02, dist02);
              if (dist02Squared <= interactionLengthSquare) break;

              const long offset1 = utils::ThreeDimensionalMapping::threeToOneD(
                  x1, y1, z1, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));

              const long offset2 = utils::ThreeDimensionalMapping::threeToOneD(
                  x2, y2, z2, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));

              // sorting direction towards middle of cell 1 and cell 2
              const std::array<double, 3> sortDirection = dist01 + dist02;

              _cellOffsets.emplace_back(offset1, offset2, utils::ArrayMath::normalize(sortDirection));
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC01Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
  ParticleCell &baseCell = cells[baseIndex];
  const size_t cOffSize = _cellOffsets.size();

  for (auto const &[offset1, offset2, r] : _cellOffsets) {
    const unsigned long otherIndex1 = baseIndex + offset1;
    const unsigned long otherIndex2 = baseIndex + offset2;
    ParticleCell &otherCell1 = cells[otherIndex1];
    ParticleCell &otherCell2 = cells[otherIndex2];

    if (baseIndex == otherIndex1 && baseIndex == otherIndex2) {
      this->_cellFunctor.processCell(baseCell);
    } else if (baseIndex == otherIndex1 && baseIndex != otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherCell2);
    } else if (baseIndex != otherIndex1 && baseIndex == otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherIndex1);
    } else {
      this->_cellFunctor.processCellTriple(baseCell, otherCell1, otherCell2);
    }
  }
}

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC01Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::resizeBuffers() {
  const auto numThreads = static_cast<size_t>(autopas_get_max_threads());
  if (_combinationSlices.size() != numThreads) {
    _combinationSlices.resize(numThreads);
    const auto cellOffsetsSize = _cellOffsets.size();
    std::for_each(_combinationSlices.begin(), _combinationSlices.end(),
                  [cellOffsetsSize](auto &e) { e.resize(cellOffsetsSize); });
    _currentSlices.resize(numThreads * _cacheOffset);
  }
}

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC01Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::traverseParticleTriplets() {
  auto &cells = *(this->_cells);
  if (not this->isApplicable()) {
    utils::ExceptionHandler::exception(
        "The C01 traversal cannot work with enabled newton3 (unless only one thread is used)!");
  }

  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
