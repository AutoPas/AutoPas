/**
 * @file RemainderTraversalEstimateManager.h
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
 * Manages the estimation of the remainder traversal time based on the number of particles in the buffer.
 */
class RemainderTraversalEstimateManager {
 public:
  /**
   * Default constructor.
   */
  RemainderTraversalEstimateManager() {
    _numParticlesInBuffer.reserve(_maxSamples);
    _remainderTraverlsalTime.reserve(_maxSamples);
  };

  /**
   * Maintains a ring buffer over distinct _numParticlesInBuffer values.
   * When the buffer is full and a new (larger) numParticleBuffer value is added, the oldest sample is removed.
   * The accumulated _sumRemainderTraversalTime is not treated as a ring buffer.
   * @param numParticleBuffer The number of particles in the buffer for which the remainder traversal time was measured.
   * @param remainderTraversalTime The measured remainder traversal time for the given numParticleBuffer.
   */
  void addNewPoint(const long numParticleBuffer, const long remainderTraversalTime) {
    _sumRemainderTraversalTime += remainderTraversalTime;

    if (_sumRemainderTraversalTime < remainderTraversalTime) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Overflow detected in remainder traversal time manager.");
    }

    // Ensures that added points are monotonically increasing in numParticleBuffer, which should be the case.
    if (_currentMaxNumParticlesInBuffer < numParticleBuffer) {
      _currentMaxNumParticlesInBuffer = numParticleBuffer;
      if (_numDifferentXConsidered < _maxSamples) {
        _numParticlesInBuffer.push_back(numParticleBuffer);
        _remainderTraverlsalTime.push_back(remainderTraversalTime);
        _numDifferentXConsidered++;
      } else {
        _numParticlesInBuffer[_currentIndex] = numParticleBuffer;
        _remainderTraverlsalTime[_currentIndex] = remainderTraversalTime;
        _currentIndex = (_currentIndex + 1) % _maxSamples;
      }
    }
  }

  /**
   * @return _sumRemainderTraversalTime, the sum of all remainder traversal times.
   */
  [[nodiscard]] long getSumRemainderTraversalTime() const { return _sumRemainderTraversalTime; }

  /**
   * Resets the manager to the empty state.
   */
  void reset() {
    _numParticlesInBuffer.clear();
    _remainderTraverlsalTime.clear();
    _currentIndex = 0;
    _numDifferentXConsidered = 0;
    _sumRemainderTraversalTime = 0;
    _currentMaxNumParticlesInBuffer = LONG_MIN;
  }

  /**
   * Predict the next remainder traversal time based on the current samples of the number of particles in the buffer and
   * their corresponding remainder traversal times.
   * @param numParticlesInBuffer The number of particles in the buffer for which the remainder traversal time should be
   * predicted.
   * @return The predicted remainder traversal time for the given numParticlesInBuffer.
   */
  [[nodiscard]] double predict(const long numParticlesInBuffer) const {
    if (_numDifferentXConsidered < _minSamples) {
      return 0;
    }
    std::pair<double, double> linRegCoeffs;
    linRegCoeffs = utils::Math::simpleLeastSquares(_numParticlesInBuffer, _remainderTraverlsalTime);
    auto prediction = linRegCoeffs.first + linRegCoeffs.second * static_cast<double>(numParticlesInBuffer);
    return prediction < 0 ? 0 : prediction;  // Remainder traversal time cannot be negative.
  }

 private:
  /**
   * Vector storing the current samples of the number of particles in the buffer for past iterations.
   * These are used for the linear regression to predict the remainder traversal time in the next iteration.
   */
  std::vector<long> _numParticlesInBuffer;

  /**
   * Vector storing the current samples of the remainder traversal time for past iterations.
   * These are used for the linear regression to predict the remainder traversal time in the next iteration.
   */
  std::vector<long> _remainderTraverlsalTime;

  /**
   * Stores the current index where the next sample will be added in the _numParticlesInBuffer and
   * _remainderTraverlsalTime vectors.
   */
  size_t _currentIndex{0};

  /**
   * Minimum number of samples needed to make a prediction.
   * Set to 2 as a linear regression requires at least two points to define a line.
   */
  const size_t _minSamples{2};

  /**
   * Maximum number of samples considered to make a prediction.
   * Set to 5 based on trials with different numbers of samples to ensure that the estimate is reliable while keeping
   * computation time low.
   */
  const size_t _maxSamples{5};

  /**
   * Counter to keep track of different number of buffer particles considered for prediction.
   */
  size_t _numDifferentXConsidered{0};

  /**
   * To keep a running sum of the remainder traversal times for all the iterations, including where the particle buffer
   * size is zero.
   */
  long _sumRemainderTraversalTime{0};

  /**
   * Used to track the largest _numParticlesInBuffer seen so far.
   * To ensure that incoming _numParticlesInBuffer values are non-decreasing.
   */
  long _currentMaxNumParticlesInBuffer = LONG_MIN;
};
}  // namespace autopas::utils::DynamicRebuildDecisionWithFPB