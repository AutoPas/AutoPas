/**
 * @file SlicedTraversal.h
 *
 * @date 20 Apr 2018
 * @author gratl
 */

#pragma once

#include <algorithm>
#include "autopas/containers/cellPairTraversals/C08BasedTraversal.h"
#include "autopas/containers/cellPairTraversals/VerletListsCellsTraversal.h"
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
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class SlicedTraversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>,
                        public VerletListsCellsTraversal<PairwiseFunctor, useNewton3> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit SlicedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>(dims, pairwiseFunctor),
        VerletListsCellsTraversal<PairwiseFunctor, useNewton3>(pairwiseFunctor) {
    rebuild(dims);
  }
  // documentation in base class
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  /**
   * Traverse verlet lists of all cells.
   * This function needs to be implemented by derived classes.
   * @param verlet verlet lists for each cell
   */
  template <class Particle>
  void traverseCellVerlet(std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>> &verlet);
  TraversalOptions getTraversalType() override;
  bool isApplicable() override;
  void rebuild(const std::array<unsigned long, 3> &dims) override;

 private:
  /**
   * Store ids of dimensions ordered by number of cells per dimensions.
   */
  std::array<int, 3> _dimsPerLength;

  /**
   * The number of cells per slice in the dimension that was sliced.
   */
  std::vector<unsigned long> _sliceThickness;
  std::vector<autopas_lock_t *> locks;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline TraversalOptions SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::getTraversalType() {
  return TraversalOptions::sliced;
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline bool SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::isApplicable() {
  return this->_sliceThickness.size() > 0;
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::rebuild(
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

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  using std::array;

  auto numSlices = _sliceThickness.size();
  // 0) check if applicable

  for (size_t i = 0; i < numSlices - 1; ++i) {
    locks[i] = new autopas_lock_t;
    autopas_init_lock(locks[i]);
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

    // all but the first slice need to lock their starting layer.
    if (slice > 0) {
      autopas_set_lock(locks[slice - 1]);
    }
    for (unsigned long dimSlice = myStartArray[_dimsPerLength[0]];
         dimSlice < myStartArray[_dimsPerLength[0]] + _sliceThickness[slice]; ++dimSlice) {
      // at the last layer request lock for the starting layer of the next
      // slice. Does not apply for the last slice.
      if (slice != numSlices - 1 && dimSlice == myStartArray[_dimsPerLength[0]] + _sliceThickness[slice] - 1) {
        autopas_set_lock(locks[slice]);
      }
      for (unsigned long dimMedium = 0; dimMedium < this->_cellsPerDimension[_dimsPerLength[1]] - 1; ++dimMedium) {
        for (unsigned long dimShort = 0; dimShort < this->_cellsPerDimension[_dimsPerLength[2]] - 1; ++dimShort) {
          array<unsigned long, 3> idArray;
          idArray[_dimsPerLength[0]] = dimSlice;
          idArray[_dimsPerLength[1]] = dimMedium;
          idArray[_dimsPerLength[2]] = dimShort;
          auto id = utils::ThreeDimensionalMapping::threeToOneD(idArray, this->_cellsPerDimension);
          this->processBaseCell(cells, id);
        }
      }
      // at the end of the first layer release the lock
      if (slice > 0 && dimSlice == myStartArray[_dimsPerLength[0]]) {
        autopas_unset_lock(locks[slice - 1]);
      } else if (slice != numSlices - 1 && dimSlice == myStartArray[_dimsPerLength[0]] + _sliceThickness[slice] - 1) {
        autopas_unset_lock(locks[slice]);
      }
    }
  }

  for (size_t i = 0; i < numSlices - 1; ++i) {
    autopas_destroy_lock(locks[i]);
    delete locks[i];
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
template <class Particle>
inline void SlicedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellVerlet(
    std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>> &verlet) {
  using std::array;

  auto numSlices = _sliceThickness.size();
  // 0) check if applicable

  for (size_t i = 0; i < numSlices - 1; ++i) {
    locks[i] = new autopas_lock_t;
    autopas_init_lock(locks[i]);
  }

#ifdef AUTOPAS_OPENMP
// although every thread gets exactly one iteration (=slice) this is faster than
// a normal parallel region
#pragma omp parallel for schedule(static, 1)
#endif
  for (size_t slice = 0; slice < numSlices; ++slice) {
    array<unsigned long, 3> myStartArray{0, 0, 0};
    for (size_t i = 0; i < slice; ++i) {
      myStartArray[_dimsPerLength[0]] += _sliceThickness[i];
    }

    // all but the first slice need to lock their starting layer.
    if (slice > 0) {
      autopas_set_lock(locks[slice - 1]);
    }
    for (unsigned long dimSlice = myStartArray[_dimsPerLength[0]];
         dimSlice < myStartArray[_dimsPerLength[0]] + _sliceThickness[slice]; ++dimSlice) {
      // at the last layer request lock for the starting layer of the next
      // slice. Does not apply for the last slice.
      if (slice != numSlices - 1 && dimSlice == myStartArray[_dimsPerLength[0]] + _sliceThickness[slice] - 1) {
        autopas_set_lock(locks[slice]);
      }
      for (unsigned long dimMedium = 0; dimMedium < this->_cellsPerDimension[_dimsPerLength[1]] - 1; ++dimMedium) {
        for (unsigned long dimShort = 0; dimShort < this->_cellsPerDimension[_dimsPerLength[2]] - 1; ++dimShort) {
          array<unsigned long, 3> idArray;
          idArray[_dimsPerLength[0]] = dimSlice;
          idArray[_dimsPerLength[1]] = dimMedium;
          idArray[_dimsPerLength[2]] = dimShort;
          auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(idArray, this->_cellsPerDimension);
          this->iterateVerletListsCell(verlet, baseIndex);
        }
      }
      // at the end of the first layer release the lock
      if (slice > 0 && dimSlice == myStartArray[_dimsPerLength[0]]) {
        autopas_unset_lock(locks[slice - 1]);
      } else if (slice != numSlices - 1 && dimSlice == myStartArray[_dimsPerLength[0]] + _sliceThickness[slice] - 1) {
        autopas_unset_lock(locks[slice]);
      }
    }
  }

  for (size_t i = 0; i < numSlices - 1; ++i) {
    autopas_destroy_lock(locks[i]);
    delete locks[i];
  }
}

}  // namespace autopas
