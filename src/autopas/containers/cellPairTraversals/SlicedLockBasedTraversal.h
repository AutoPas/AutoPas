/**
 * @file SlicedLockBasedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include "autopas/containers/cellPairTraversals/SlicedBasedTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into multiple slices along this dimension. Slices are
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
class SlicedLockBasedTraversal : public SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit SlicedLockBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                    const double interactionLength, const std::array<double, 3> &cellLength)
      : SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                    interactionLength, cellLength) {}

 protected:
  /**
   * whether to use static or dynamic scheduling.
   */
  bool _dynamic = true;

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
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <bool allCells, typename LoopBody>
void SlicedLockBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::slicedTraversal(
    LoopBody &&loopBody) {
  using std::array;

  auto numSlices = this->_sliceThickness.size();
  std::vector<AutoPasLock> locks;
  locks.resize((numSlices - 1) * this->_overlapLongestAxis);

  // 0) check if applicable
  std::array<size_t, 2> overLapps23{this->_overlap[this->_dimsPerLength[1]], this->_overlap[this->_dimsPerLength[2]]};
  if (allCells) {
    overLapps23 = {0ul, 0ul};
    this->_sliceThickness.back() += this->_overlapLongestAxis;
  }

  std::vector<utils::Timer> timers;
  std::vector<double> threadTimes;

  timers.resize(numSlices);
  threadTimes.resize(numSlices);

#ifdef AUTOPAS_OPENMP
  // although every thread gets exactly one iteration (=slice) this is faster than a normal parallel region
  auto numThreads = static_cast<size_t>(autopas_get_max_threads());
  if (this->_dynamic) {
    omp_set_schedule(omp_sched_dynamic, 1);
  } else {
    omp_set_schedule(omp_sched_static, 1);
  }
#pragma omp parallel for schedule(runtime) num_threads(numThreads)
#endif
  for (size_t slice = 0; slice < numSlices; ++slice) {
    timers[slice].start();
    array<unsigned long, 3> myStartArray{0, 0, 0};
    for (size_t i = 0; i < slice; ++i) {
      myStartArray[this->_dimsPerLength[0]] += this->_sliceThickness[i];
    }

    // all but the first slice need to lock their starting layers.
    const unsigned long lockBaseIndex = (slice - 1) * this->_overlapLongestAxis;
    if (slice > 0) {
      for (unsigned long i = 0ul; i < this->_overlapLongestAxis; i++) {
        locks[lockBaseIndex + i].lock();
      }
    }
    const auto lastLayer = myStartArray[this->_dimsPerLength[0]] + this->_sliceThickness[slice];
    for (unsigned long sliceOffset = 0ul; sliceOffset < this->_sliceThickness[slice]; ++sliceOffset) {
      const auto dimSlice = myStartArray[this->_dimsPerLength[0]] + sliceOffset;
      // at the last layers request lock for the starting layer of the next
      // slice. Does not apply for the last slice.
      if (slice != numSlices - 1 and dimSlice >= lastLayer - this->_overlapLongestAxis) {
        locks[((slice + 1) * this->_overlapLongestAxis) - (lastLayer - dimSlice)].lock();
      }
      for (unsigned long dimMedium = 0; dimMedium < this->_cellsPerDimension[this->_dimsPerLength[1]] - overLapps23[0];
           ++dimMedium) {
        for (unsigned long dimShort = 0; dimShort < this->_cellsPerDimension[this->_dimsPerLength[2]] - overLapps23[1];
             ++dimShort) {
          array<unsigned long, 3> idArray = {};
          idArray[this->_dimsPerLength[0]] = dimSlice;
          idArray[this->_dimsPerLength[1]] = dimMedium;
          idArray[this->_dimsPerLength[2]] = dimShort;

          loopBody(idArray[0], idArray[1], idArray[2]);
        }
      }
      // at the end of the first layers release the lock
      if (slice > 0 and dimSlice < myStartArray[this->_dimsPerLength[0]] + this->_overlapLongestAxis) {
        locks[lockBaseIndex + sliceOffset].unlock();
        // if lastLayer is reached within overlap area, unlock all following locks
        // this should never be the case if slice thicknesses are set up properly; thickness should always be
        // greater than the overlap along the longest axis, or the slices won't be processed in parallel.
        if (dimSlice == lastLayer - 1) {
          for (unsigned long i = sliceOffset + 1; i < this->_overlapLongestAxis; ++i) {
            locks[lockBaseIndex + i].unlock();
          }
        }
      } else if (slice != numSlices - 1 and dimSlice == lastLayer - 1) {
        // clearing of the locks set on the last layers of each slice
        for (size_t i = (slice * this->_overlapLongestAxis); i < (slice + 1) * this->_overlapLongestAxis; ++i) {
          locks[i].unlock();
        }
      }
    }
    threadTimes[slice] = timers[slice].stop();
  }

  if (allCells) {
    this->_sliceThickness.back() -= this->_overlapLongestAxis;
  }

  std::string timesStr;
  for (auto t : threadTimes) {
    timesStr += std::to_string(t) + ", ";
  }
  auto minMax = std::minmax_element(threadTimes.begin(), threadTimes.end());
  auto avg = std::accumulate(threadTimes.begin(), threadTimes.end(), 0.0) / numSlices;
  auto variance = std::accumulate(threadTimes.cbegin(), threadTimes.cend(), 0.0,
                                  [avg](double a, double b) -> double { return a + std::pow(avg - b, 2.0); }) /
                  numSlices;
  auto stddev = std::sqrt(variance);

  AutoPasLog(debug, "times per slice: [{}].", timesStr);
  AutoPasLog(debug, "Difference between longest and shortest time: {:.3G}", *minMax.second - *minMax.first);
  AutoPasLog(debug, "Ratio between longest and shortest time: {:.3G}", (float)*minMax.second / *minMax.first);
  AutoPasLog(debug, "avg: {:.3G}, std-deviation: {:.3G} ({:.3G}%)", avg, stddev, 100 * stddev / avg);
}

}  // namespace autopas
