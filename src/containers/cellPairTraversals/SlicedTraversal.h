/*
 * SlicedTraversal.h
 *
 *  Created on: 20 Apr 2018
 *      Author: gratl
 */

#pragma once

#include <utils/WrapOpenMP.h>
#include <algorithm>
#include "C08BasedTraversal.h"
#include "utils/ThreeDimensionalMapping.h"

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
 * @tparam ParticleCell the type of cells
 * @tparam CellFunctor the cell functor that defines the interaction of the
 * particles of two specific cells
 */
template <class ParticleCell, class CellFunctor>
class SlicedTraversal : public C08BasedTraversal<ParticleCell, CellFunctor> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param cells the cells through which the traversal should traverse.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param cellfunctor The cell functor that defines the interaction of
   * particles between two different cells.
   */
  explicit SlicedTraversal(std::vector<ParticleCell> &cells,
                           const std::array<unsigned long, 3> &dims,
                           CellFunctor *cellfunctor)
      : C08BasedTraversal<ParticleCell, CellFunctor>(cells, dims, cellfunctor) {
    rebuild(cells, dims);
  }
  // documentation in base class
  void traverseCellPairs() override;
  bool isApplicable() override;
  void rebuild(std::vector<ParticleCell> &cells,
               const std::array<unsigned long, 3> &dims) override;

 private:
  // FIXME: Remove this as soon as other traversals are available
  void traverseCellPairsFallback();

  /**
   * store ids of dimensions ordered by number of cells per dimensions
   */
  std::array<int, 3> _dimsPerLength;

  /**
   * the number of cells per slice in the dimension that was slicedFjjkj
   */
  std::vector<unsigned long> _sliceThickness;
  std::vector<autopas_lock_t *> locks;
};

template <class ParticleCell, class CellFunctor>
inline bool SlicedTraversal<ParticleCell, CellFunctor>::isApplicable() {
  return this->_cellsPerDimension[_dimsPerLength[0]] /
             this->_sliceThickness.size() >=
         2;
}

template <class ParticleCell, class CellFunctor>
inline void SlicedTraversal<ParticleCell, CellFunctor>::rebuild(
    std::vector<ParticleCell> &cells,
    const std::array<unsigned long, 3> &dims) {
  CellPairTraversals<ParticleCell, CellFunctor>::rebuild(cells, dims);

  // find longest dimension
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(),
                                        this->_cellsPerDimension.end());
  _dimsPerLength[0] =
      (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsPerLength[2] =
      (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsPerLength[1] = 3 - (_dimsPerLength[0] + _dimsPerLength[2]);

  // split domain across its longest dimension

  auto numSlices = autopas_get_max_threads();
  // find dimension with most cells
  _sliceThickness.clear();
  _sliceThickness.insert(
      _sliceThickness.begin(), numSlices,
      this->_cellsPerDimension[_dimsPerLength[0]] / numSlices);
  auto rest = this->_cellsPerDimension[_dimsPerLength[0]] -
              _sliceThickness[0] * numSlices;
  for (int i = 0; i < rest; ++i) ++_sliceThickness[i];
  // decreases last _sliceThickness by one to account for the way we handle base
  // cells
  --*--_sliceThickness.end();

  locks.resize(numSlices);
}

// FIXME: Remove this as soon as other traversals are available
template <class ParticleCell, class CellFunctor>
inline void
SlicedTraversal<ParticleCell, CellFunctor>::traverseCellPairsFallback() {
  std::array<unsigned long, 3> endid;
  for (int d = 0; d < 3; ++d) {
    endid[d] = this->_cellsPerDimension[d] - 1;
  }
  for (unsigned long z = 0; z < endid[2]; ++z) {
    for (unsigned long y = 0; y < endid[1]; ++y) {
      for (unsigned long x = 0; x < endid[0]; ++x) {
        unsigned long ind = ThreeDimensionalMapping::threeToOneD(
            x, y, z, this->_cellsPerDimension);
        this->processBaseCell(ind);
      }
    }
  }
}

template <class ParticleCell, class CellFunctor>
inline void SlicedTraversal<ParticleCell, CellFunctor>::traverseCellPairs() {
  using std::array;

  auto numSlices = _sliceThickness.size();
  // 0) check if applicable

  // use fallback version if not applicable
  // FIXME: Remove this as soon as other traversals are available
  if (not isApplicable()) {
    traverseCellPairsFallback();
    return;
  }

  for (auto i = 0; i < numSlices - 1; ++i) {
    locks[i] = new autopas_lock_t;
    autopas_init_lock(locks[i]);
  }

#ifdef AUTOPAS_OPENMP
// although every thread gets exactly one iteration (=slice) this is faster than
// a normal parallel region
#pragma omp parallel for schedule(static, 1)
#endif
  for (auto slice = 0; slice < numSlices; ++slice) {
    array<unsigned long, 3> myStartArray{0, 0, 0};
    for (int i = 0; i < slice; ++i) {
      myStartArray[_dimsPerLength[0]] += _sliceThickness[i];
    }

    // all but the first slice need to lock their starting layer.
    if (slice > 0) {
      autopas_set_lock(locks[slice - 1]);
    }
    for (unsigned long dimSlice = myStartArray[_dimsPerLength[0]];
         dimSlice < myStartArray[_dimsPerLength[0]] + _sliceThickness[slice];
         ++dimSlice) {
      // at the last layer request lock for the starting layer of the next
      // slice. Does not apply for the last slice.
      if (slice != numSlices - 1 &&
          dimSlice ==
              myStartArray[_dimsPerLength[0]] + _sliceThickness[slice] - 1) {
        autopas_set_lock(locks[slice]);
      }
      for (unsigned long dimMedium = 0;
           dimMedium < this->_cellsPerDimension[_dimsPerLength[1]] - 1;
           ++dimMedium) {
        for (unsigned long dimShort = 0;
             dimShort < this->_cellsPerDimension[_dimsPerLength[2]] - 1;
             ++dimShort) {
          array<unsigned long, 3> idArray;
          idArray[_dimsPerLength[0]] = dimSlice;
          idArray[_dimsPerLength[1]] = dimMedium;
          idArray[_dimsPerLength[2]] = dimShort;
          auto id = ThreeDimensionalMapping::threeToOneD(
              idArray, this->_cellsPerDimension);
          this->processBaseCell(id);
        }
      }
      // at the end of the first layer release the lock
      if (slice > 0 && dimSlice == myStartArray[_dimsPerLength[0]]) {
        autopas_unset_lock(locks[slice - 1]);
      } else if (slice != numSlices - 1 &&
                 dimSlice == myStartArray[_dimsPerLength[0]] +
                                 _sliceThickness[slice] - 1) {
        autopas_unset_lock(locks[slice]);
      }
    }
  }

  for (auto i = 0; i < numSlices - 1; ++i) {
    autopas_destroy_lock(locks[i]);
    delete locks[i];
  }
}

} /* namespace autopas */