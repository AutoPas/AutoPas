/**
 * @file SlicedLockBasedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include <numeric>

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
 * @tparam spaciallyForward Whether the base step only covers neigboring cells tha are spacially forward (for example
 * c08)
 */
template <class ParticleCell, class PairwiseFunctor, bool spaciallyForward>
class SlicedLockBasedTraversal : public SlicedBasedTraversal<ParticleCell, PairwiseFunctor, spaciallyForward> {
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
                                    const double interactionLength, const std::array<double, 3> &cellLength,
                                    DataLayoutOption::Value dataLayout, bool useNewton3)
      : SlicedBasedTraversal<ParticleCell, PairwiseFunctor, spaciallyForward>(dims, pairwiseFunctor, interactionLength,
                                                                              cellLength, dataLayout, useNewton3) {}

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
   */
  template <typename LoopBody>
  inline void slicedTraversal(LoopBody &&loopBody);
};

template <class ParticleCell, class PairwiseFunctor, bool spaciallyForward>
template <typename LoopBody>
void SlicedLockBasedTraversal<ParticleCell, PairwiseFunctor, spaciallyForward>::slicedTraversal(LoopBody &&loopBody) {
  using std::array;

  auto numSlices = this->_sliceThickness.size();
  std::vector<AutoPasLock> locks;
  locks.resize((numSlices - 1) * this->_overlapLongestAxis);

  // 0) check if applicable
  const auto overLapps23 = [&]() -> std::array<size_t, 2> {
    if constexpr (spaciallyForward) {
      return {this->_overlap[this->_dimsPerLength[1]], this->_overlap[this->_dimsPerLength[2]]};
    } else {
      return {0ul, 0ul};
    }
  }();

  std::vector<utils::Timer> timers;
  std::vector<double> threadTimes;

  timers.resize(numSlices);
  threadTimes.resize(numSlices);

#ifdef AUTOPAS_USE_OPENMP
  if (this->_dynamic) {
    omp_set_schedule(omp_sched_dynamic, 1);
  } else {
    omp_set_schedule(omp_sched_static, 1);
  }
#endif
  AUTOPAS_OPENMP(parallel for schedule(runtime))
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

  AutoPasLog(DEBUG, "times per slice: [{}].", timesStr);
  AutoPasLog(DEBUG, "Difference between longest and shortest time: {:.3G}", *minMax.second - *minMax.first);
  AutoPasLog(DEBUG, "Ratio between longest and shortest time: {:.3G}", (float)*minMax.second / *minMax.first);
  AutoPasLog(DEBUG, "avg: {:.3G}, std-deviation: {:.3G} ({:.3G}%)", avg, stddev, 100 * stddev / avg);
}

}  // namespace autopas
