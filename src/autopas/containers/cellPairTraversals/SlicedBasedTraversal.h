/**
 * @file SlicedBasedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain in one slice (block) per thread along this dimension. Slices are
 * assigned to the threads in a round robin fashion. Each thread locks the cells
 * on the boundary wall to the previous slice with one lock. This lock is lifted
 * as soon the boundary wall is fully processed.
 *
 * @tparam ParticleCell The type of cells.
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class SlicedBasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit SlicedBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                const double interactionLength, const std::array<double, 3> &cellLength)
      : CellPairTraversal<ParticleCell>(dims),
        _overlap{},
        _dimsPerLength{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlapLongestAxis(0),
        _sliceThickness{},
        locks(),
        _dataLayoutConverter(pairwiseFunctor) {
    init(dims);
  }

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  bool isApplicable() const override {
    return not(dataLayout == DataLayoutOption::cuda) and this->_sliceThickness.size() > 0;
  }

  /**
   * Load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  void initTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      /// @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.loadDataLayout(cells[i]);
      }
    }
  }

  /**
   * Write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      /// @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.storeDataLayout(cells[i]);
      }
    }
  }

 protected:
  /**
   * Resets the cell structure of the traversal.
   * @param dims
   */
  void init(const std::array<unsigned long, 3> &dims);

  /**
   * The main traversal of the sliced traversal.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   *
   * @tparam allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the sliced step if allCells is false, iteration will not occur over the last layer of
   * cells (for _overlap=1) (in x, y and z direction).
   */
  template <bool allCells = false, typename LoopBody>
  inline void slicedTraversal(LoopBody &&loopBody);

  /**
   * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  std::array<unsigned long, 3> _overlap;

  /**
   * Interaction length (cutoff + skin).
   */
  double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  std::array<double, 3> _cellLength;

 private:
  /**
   * Store ids of dimensions ordered by number of cells per dimensions.
   */
  std::array<int, 3> _dimsPerLength;

  /**
   * Overlap of interacting cells along the longest axis.
   */
  unsigned long _overlapLongestAxis;

  /**
   * The number of cells per slice in the dimension that was sliced.
   */
  std::vector<unsigned long> _sliceThickness;
  std::vector<AutoPasLock> locks;

  /**
   * Data Layout Converter to be used with this traversal.
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::init(
    const std::array<unsigned long, 3> &dims) {
  for (unsigned int d = 0; d < 3; d++) {
    _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
  }

  // find longest dimension
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(), this->_cellsPerDimension.end());
  _dimsPerLength[0] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsPerLength[2] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsPerLength[1] = 3 - (_dimsPerLength[0] + _dimsPerLength[2]);

  _overlapLongestAxis = _overlap[_dimsPerLength[0]];

  // split domain across its longest dimension

  auto numSlices = (size_t)autopas_get_max_threads();
  auto minSliceThickness = this->_cellsPerDimension[_dimsPerLength[0]] / numSlices;
  if (minSliceThickness < _overlapLongestAxis + 1) {
    minSliceThickness = _overlapLongestAxis + 1;
    numSlices = this->_cellsPerDimension[_dimsPerLength[0]] / minSliceThickness;
    AutoPasLog(debug, "Sliced traversal only using {} threads because the number of cells is too small.", numSlices);
  }

  _sliceThickness.clear();

  // abort if domain is too small -> cleared _sliceThickness array indicates non applicability
  if (numSlices < 1) return;

  _sliceThickness.insert(_sliceThickness.begin(), numSlices, minSliceThickness);
  auto rest = this->_cellsPerDimension[_dimsPerLength[0]] - _sliceThickness[0] * numSlices;
  for (size_t i = 0; i < rest; ++i) ++_sliceThickness[i];
  // decreases last _sliceThickness by _overlapLongestAxis to account for the way we handle base cells
  _sliceThickness.back() -= _overlapLongestAxis;

  locks.resize((numSlices - 1) * _overlapLongestAxis);
}
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <bool allCells, typename LoopBody>
void SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::slicedTraversal(LoopBody &&loopBody) {
  using std::array;

  auto numSlices = _sliceThickness.size();
  // 0) check if applicable

  std::array<size_t, 2> overLapps23{_overlap[_dimsPerLength[1]], _overlap[_dimsPerLength[2]]};
  if (allCells) {
    overLapps23 = {0ul, 0ul};
    _sliceThickness.back() += _overlapLongestAxis;
  }

#ifdef AUTOPAS_OPENMP
// although every thread gets exactly one iteration (=slice) this is faster than a normal parallel region
#pragma omp parallel for schedule(static, 1) num_threads(numSlices)
#endif
  for (size_t slice = 0; slice < numSlices; ++slice) {
    array<unsigned long, 3> myStartArray{0, 0, 0};
    for (size_t i = 0; i < slice; ++i) {
      myStartArray[_dimsPerLength[0]] += _sliceThickness[i];
    }

    // all but the first slice need to lock their starting layers.
    if (slice > 0) {
      for (unsigned long i = 0ul; i < _overlapLongestAxis; i++) {
        locks[((slice - 1) * _overlapLongestAxis) + i].lock();
      }
    }
    const auto lastLayer = myStartArray[_dimsPerLength[0]] + _sliceThickness[slice];
    for (unsigned long dimSlice = myStartArray[_dimsPerLength[0]]; dimSlice < lastLayer; ++dimSlice) {
      // at the last layers request lock for the starting layer of the next
      // slice. Does not apply for the last slice.
      if (slice != numSlices - 1 && dimSlice >= lastLayer - _overlapLongestAxis) {
        locks[((slice + 1) * _overlapLongestAxis) - (lastLayer - dimSlice)].lock();
      }
      for (unsigned long dimMedium = 0; dimMedium < this->_cellsPerDimension[_dimsPerLength[1]] - overLapps23[0];
           ++dimMedium) {
        for (unsigned long dimShort = 0; dimShort < this->_cellsPerDimension[_dimsPerLength[2]] - overLapps23[1];
             ++dimShort) {
          array<unsigned long, 3> idArray = {};
          idArray[_dimsPerLength[0]] = dimSlice;
          idArray[_dimsPerLength[1]] = dimMedium;
          idArray[_dimsPerLength[2]] = dimShort;
          loopBody(idArray[0], idArray[1], idArray[2]);
        }
      }
      // at the end of the first layers release the lock
      if (slice > 0 && dimSlice < myStartArray[_dimsPerLength[0]] + _overlapLongestAxis) {
        const unsigned long index = ((slice - 1) * _overlapLongestAxis) + (dimSlice - myStartArray[_dimsPerLength[0]]);
        locks[index].unlock();
        // if lastLayer is reached within overlap area, unlock all following locks
        if (dimSlice == lastLayer - 1) {
          for (unsigned long i = dimSlice + 1; i < myStartArray[_dimsPerLength[0]] + _overlapLongestAxis; ++i) {
            const unsigned long index = ((slice - 1) * _overlapLongestAxis) + (i - myStartArray[_dimsPerLength[0]]);
            locks[index].unlock();
          }
        }
      } else if (slice != numSlices - 1 && dimSlice == lastLayer - 1) {
        // clearing of the locks set on the last layers of each slice
        for (size_t i = (slice * _overlapLongestAxis); i < (slice + 1) * _overlapLongestAxis; ++i) {
          locks[i].unlock();
        }
      }
    }
  }

  if (allCells) {
    _sliceThickness.back() -= _overlapLongestAxis;
  }
}

}  // namespace autopas
