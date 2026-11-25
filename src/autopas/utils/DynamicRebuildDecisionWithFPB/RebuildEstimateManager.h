/**
 * @file RebuildEstimateManager.h
 * @author Manish Mishra
 * @date 25.11.2026
 */

#pragma once

#include <algorithm>
#include <climits>
#include <iostream>
#include <limits>
#include <ostream>
#include <utility>
#include <vector>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"

namespace autopas::utils::DynamicRebuildDecisionWithFPB {

/**
 * @brief Base class for managing estimation of rebuild time for dynamic rebuilding with fast particles buffer.
 * The class maintains a ring buffer of the last five rebuild time samples and provides methods to add new samples and
 * predict the rebuild time based on the current samples.
 */
class RebuildEstimateManager {
 public:
  /**
   * Default constructor.
   */
  RebuildEstimateManager() { _rebuildTimes.assign(_numMaxSamples, 0); };

  /**
   * Adds a new rebuild time sample to the manager.
   * @param rebuildTime The rebuild time sample to add.
   */
  void addNewPoint(long const rebuildTime) {
    _rebuildTimesSum -= _rebuildTimes[_currentIndex];
    _rebuildTimes[_currentIndex] = rebuildTime;
    _rebuildTimesSum += rebuildTime;

    if (_rebuildTimesSum < rebuildTime) {
      throw autopas::utils::ExceptionHandler::AutoPasException("Overflow detected in rebuild time manager.");
    }

    _currentIndex = (_currentIndex + 1) % _numMaxSamples;
    _numSamples++;
  }

  /**
   * Predicts the rebuild time based on the current samples.
   * @return The predicted rebuild time.
   */
  double predict() {
    if (_numSamples == 0) {
      throw autopas::utils::ExceptionHandler::AutoPasException("No rebuild samples available for prediction.");
    }
    auto rebuildEstimate = static_cast<double>(_rebuildTimesSum) / std::min(_numMaxSamples, _numSamples);
    return rebuildEstimate;
  }

  /**
   * Resets the manager to the empty state.
   */
  void reset() {
    _rebuildTimes.assign(_numMaxSamples, 0);
    _currentIndex = 0;
    _numSamples = 0;
    _rebuildTimesSum = 0;
  }

 private:
  /**
   * Stores upto last five samples of the rebuild times.
   */
  std::vector<long> _rebuildTimes;

  /**
   * Stores the sum of the current samples of the rebuild times.
   */
  long _rebuildTimesSum{0};

  /**
   * Stores the current number of samples in the _rebuildTimes vector.
   */
  size_t _numSamples{0};

  /**
   * Stores the current index in the _rebuildTimes vector where the next sample will be added.
   */
  size_t _currentIndex{0};

  /**
   * Stores the maximum number of samples in the _rebuildTimes vector.
   * The values is set to 5 to ensure that the estimate is reliable and reflects the current state of the simulation and
   * is not overly influenced by older samples. The magic number 5 was chosen by trial and error with different number
   * of samples.
   */
  const size_t _numMaxSamples{5};
};
}  // namespace autopas::utils::DynamicRebuildDecisionWithFPB