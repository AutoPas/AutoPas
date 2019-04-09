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
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class SlicedBasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit SlicedBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims),
        _dimsPerLength{},
        _sliceThickness{},
        locks(),
        _dataLayoutConverter(pairwiseFunctor) {
    rebuild(dims);
  }

  bool isApplicable() override {
#if defined(AUTOPAS_CUDA)
    if (DataLayout == DataLayoutOption::cuda) {
      int nDevices;
      cudaGetDeviceCount(&nDevices);
      return (this->_sliceThickness.size() > 0) && (nDevices > 0);
    } else {
      return this->_sliceThickness.size() > 0;
    }
#else
    return this->_sliceThickness.size() > 0;
#endif
  }

  void initTraversal(std::vector<ParticleCell> &cells) override {
    for (auto &cell : cells) {
      _dataLayoutConverter.loadDataLayout(cell);
    }
  }

  void endTraversal(std::vector<ParticleCell> &cells) override {
    for (auto &cell : cells) {
      _dataLayoutConverter.storeDataLayout(cell);
    }
  }
  void rebuild(const std::array<unsigned long, 3> &dims) override;

 protected:
  /**
   * The main traversal of the C01Traversal.
   * @copydetails C01BasedTraversal::c01Traversal()
   */
  template <typename LoopBody>
  inline void slicedTraversal(LoopBody &&loopBody);

 private:
  /**
   * Store ids of dimensions ordered by number of cells per dimensions.
   */
  std::array<int, 3> _dimsPerLength;

  /**
   * The number of cells per slice in the dimension that was sliced.
   */
  std::vector<unsigned long> _sliceThickness;
  std::vector<AutoPasLock> locks;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, DataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void SlicedBasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::rebuild(
    const std::array<unsigned long, 3> &dims) {
  CellPairTraversal<ParticleCell>::rebuild(dims);

  // find longest dimension
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(), this->_cellsPerDimension.end());
  _dimsPerLength[0] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsPerLength[2] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsPerLength[1] = 3 - (_dimsPerLength[0] + _dimsPerLength[2]);

  // split domain across its longest dimension

  auto numSlices = (size_t)autopas_get_max_threads();
  auto minSliceThickness = this->_cellsPerDimension[_dimsPerLength[0]] / numSlices;
  if (minSliceThickness < 2) {
    minSliceThickness = 2;
    numSlices = this->_cellsPerDimension[_dimsPerLength[0]] / minSliceThickness;
    AutoPasLog(debug, "Sliced traversal only using {} threads because the number of cells is too small.", numSlices);
  }

  _sliceThickness.clear();

  // abort if domain is too small -> cleared _sliceThickness array indicates non applicability
  if (numSlices < 1) return;

  _sliceThickness.insert(_sliceThickness.begin(), numSlices, minSliceThickness);
  auto rest = this->_cellsPerDimension[_dimsPerLength[0]] - _sliceThickness[0] * numSlices;
  for (size_t i = 0; i < rest; ++i) ++_sliceThickness[i];
  // decreases last _sliceThickness by one to account for the way we handle base cells
  --*--_sliceThickness.end();

  locks.resize(numSlices);
}
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
template <typename LoopBody>
void SlicedBasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::slicedTraversal(LoopBody &&loopBody) {
  using std::array;

  auto numSlices = _sliceThickness.size();
  // 0) check if applicable

#ifdef AUTOPAS_OPENMP
// although every thread gets exactly one iteration (=slice) this is faster than a normal parallel region
#pragma omp parallel for schedule(static, 1) num_threads(numSlices)
#endif
  for (size_t slice = 0; slice < numSlices; ++slice) {
    array<unsigned long, 3> myStartArray{0, 0, 0};
    for (size_t i = 0; i < slice; ++i) {
      myStartArray[_dimsPerLength[0]] += _sliceThickness[i];
    }

    // all but the first slice need to lock their starting layer.
    if (slice > 0) {
      locks[slice - 1].lock();
    }
    for (unsigned long dimSlice = myStartArray[_dimsPerLength[0]];
         dimSlice < myStartArray[_dimsPerLength[0]] + _sliceThickness[slice]; ++dimSlice) {
      // at the last layer request lock for the starting layer of the next
      // slice. Does not apply for the last slice.
      if (slice != numSlices - 1 && dimSlice == myStartArray[_dimsPerLength[0]] + _sliceThickness[slice] - 1) {
        locks[slice].lock();
      }
      for (unsigned long dimMedium = 0; dimMedium < this->_cellsPerDimension[_dimsPerLength[1]] - 1; ++dimMedium) {
        for (unsigned long dimShort = 0; dimShort < this->_cellsPerDimension[_dimsPerLength[2]] - 1; ++dimShort) {
          array<unsigned long, 3> idArray = {};
          idArray[_dimsPerLength[0]] = dimSlice;
          idArray[_dimsPerLength[1]] = dimMedium;
          idArray[_dimsPerLength[2]] = dimShort;
          loopBody(idArray[0], idArray[1], idArray[2]);
        }
      }
      // at the end of the first layer release the lock
      if (slice > 0 && dimSlice == myStartArray[_dimsPerLength[0]]) {
        locks[slice - 1].unlock();
      } else if (slice != numSlices - 1 && dimSlice == myStartArray[_dimsPerLength[0]] + _sliceThickness[slice] - 1) {
        // clearing of the lock set on the last layer of each slice
        locks[slice].unlock();
      }
    }
  }
}

}  // namespace autopas
