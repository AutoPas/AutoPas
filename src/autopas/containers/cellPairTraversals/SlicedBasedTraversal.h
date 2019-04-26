/**
 * @file SlicedBasedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include <algorithm>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
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
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class SlicedBasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit SlicedBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                const double cutoff = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : CellPairTraversal<ParticleCell>(dims),
        _dimsPerLength{},
        _cutoff(cutoff),
        _cellLength(cellLength),
        _overlap(0ul),
        _sliceThickness{},
        locks() {
    rebuild(dims);
  }

  bool isApplicable() override { return this->_sliceThickness.size() > 0; }

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
   * cutoff radius.
   */
  double _cutoff;

  /**
   * cell length in CellBlock3D.
   */
  std::array<double, 3> _cellLength;

  /**
   * overlap between interacting cells along longest axis
   */
  unsigned long _overlap;

  /**
   * The number of cells per slice in the dimension that was sliced.
   */
  std::vector<unsigned long> _sliceThickness;
  std::vector<AutoPasLock> locks;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void SlicedBasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::rebuild(
    const std::array<unsigned long, 3> &dims) {
  CellPairTraversal<ParticleCell>::rebuild(dims);

  // find longest dimension
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(), this->_cellsPerDimension.end());
  _dimsPerLength[0] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsPerLength[2] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsPerLength[1] = 3 - (_dimsPerLength[0] + _dimsPerLength[2]);

  _overlap = std::ceil(_cutoff / _cellLength[_dimsPerLength[0]]);

  // split domain across its longest dimension

  auto numSlices = (size_t)autopas_get_max_threads();
  auto minSliceThickness = this->_cellsPerDimension[_dimsPerLength[0]] / numSlices;
  if (minSliceThickness < _overlap + 1) {
    minSliceThickness = _overlap + 1;
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

  locks.resize(numSlices * _overlap);
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
template <typename LoopBody>
void SlicedBasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::slicedTraversal(LoopBody &&loopBody) {
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

    // all but the first slice need to lock their starting layers.
    if (slice > 0) {
      for (unsigned long i = 1ul; i <= _overlap; i++) {
        locks[(slice * _overlap) - i].lock();
      }
    }
    const auto lastLayer = myStartArray[_dimsPerLength[0]] + _sliceThickness[slice];
    for (unsigned long dimSlice = myStartArray[_dimsPerLength[0]]; dimSlice < lastLayer; ++dimSlice) {
      // at the last layer request lock for the starting layer of the next
      // slice. Does not apply for the last slice.
      if (slice != numSlices - 1 && dimSlice >= lastLayer - _overlap) {
        locks[(slice * _overlap) + _overlap - (lastLayer - dimSlice)].lock();
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
      // at the end of the first layers release the lock
      if (slice > 0 && dimSlice < myStartArray[_dimsPerLength[0]] + _overlap) {
        locks[(slice * _overlap) - (_overlap - (dimSlice - myStartArray[_dimsPerLength[0]]))].unlock();
      } else if (slice != numSlices - 1 && dimSlice == lastLayer - 1) {
        // clearing of the locks set on the last layers of each slice
        for (size_t i = (slice * _overlap); i < (slice + 1) * _overlap; ++i) {
          locks[i].unlock();
        }
      }
    }
  }
}

}  // namespace autopas
