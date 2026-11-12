/**
 * @file RegressionModels.h
 * @author Natallia Padalinskaya
 * @date 02.12.2025
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

namespace autopas::utils::RegressionModels {

/**
 * Base class for estimators that predict a long-valued response variable y
 * based on streaming runtime measurements (e.g. remainder traversal or rebuild times).
 *
 * Provides common state such as sample counters and point limits.
 * Concrete estimation strategies (mean, linear regression) are implemented
 * in derived classes.
 */
class RegressionBase {
 public:
  // TODO: minN < maxN & minN > 0
  /**
   * @param minN Lower limit for the number of points considered for prediction.
   * @param maxN Upper limit for the number of points considered for prediction.
   */
  RegressionBase(const size_t minN, const size_t maxN) : _minN(minN), _maxN(maxN) {}

  /// Return codes for prediction status.
  enum class ReturnCode { OK_REG, NOT_ENOUGH_POINTS_REG, EXCEEDED_MAX_POINTS_REG, OVERFLOW_REG, UNKNOWN_ERROR_REG };

  /**
   * Output type for the predict method, used to capture whether an error occurred during prediction, e.g., a division
   * by zero.
   */
  struct Result {
    double _value;
    ReturnCode _returnCode;
    bool _isOk;
  };

  /// Reset to initial values.
  void reset() {
    _n = 0;
    _sumY = 0.0;
  }

  /**
   * adds y to sum an increments counter of points
   * updates - _sumY and _n
   */
  ReturnCode addNewPoint(long const y) {
    _sumY = Math::safeAdd(_sumY, y);
    if (isOverflown(_sumY)) {
      return ReturnCode::OVERFLOW_REG;
    }
    incrementN();
    return ReturnCode::OK_REG;
  }

  /// Increment counter n when adding a new point.
  void incrementN() { _n++; }

  /// Checks whether a value indicates an arithmetic overflow based on sentinel limits.
  static bool isOverflown(const long x) {
    if (x == std::numeric_limits<long>::max() || x == std::numeric_limits<long>::min()) {
      return true;
    }
    return false;
  }

  /// Checks whether at least _minN points have been gathered, the minimum number required for prediction.
  [[nodiscard]] bool hasEnoughPoints() const { return _n >= _minN; }
  /// Checks whether _maxN points have already been gathered.
  [[nodiscard]] bool reachedMaxPoints() const { return _n >= _maxN; }
  /// Checks whether more than _maxN points have already been gathered.
  [[nodiscard]] bool exceedsMaxPoints() const { return _n > _maxN; }

  /// Returns the number of points currently added.
  [[nodiscard]] size_t getN() const { return _n; }

  /// Returns the sum of all added y values.
  [[nodiscard]] long getSumY() const { return _sumY; }

 protected:
  /// Number of data points currently stored
  size_t _n = 0;

  /// Sum of all y values collected so far.
  long _sumY = 0;

  /// Lower limit for the number of points considered for prediction.
  size_t _minN = 0;

  /// Upper limit for the number of points considered for prediction.
  size_t _maxN = 0;
};

/**
 * Estimates a constant value by computing the arithmetic mean of observed
 * y values in a streaming fashion.
 *
 * Used as a lightweight baseline estimator for quantities with low variance
 * (e.g. rebuildNeighborListTime).
 */
class RunningMean : public RegressionBase {
 public:
  /**
   * Default constructor.
   * Initializes the base RegressionBase with minN = 1 and maxN = 100.
   */
  RunningMean() : RegressionBase(1, std::numeric_limits<long>::max()) {}

  /// Reset to initial values.
  void reset() { RegressionBase::reset(); }

  /**
   *If the maximum number of data points has been reached,no further samples are added to keep the estimator bounded.
   */
  ReturnCode addNewPoint(long const y) { return RegressionBase::addNewPoint(y); }

  /// The prediction of the Mean class takes the mean over the gathered points.
  Result predict() {
    if (not hasEnoughPoints()) {
      return Result{0, ReturnCode::NOT_ENOUGH_POINTS_REG, false};
    }

    /*
     * If the maximum number of points has not been reached yet,
     * recompute the mean using the current accumulated values.
     * Once the estimator is full, the cached mean can be reused.
     */
    auto runningMean = static_cast<double>(_sumY) / static_cast<double>(_n);
    return Result{runningMean, ReturnCode::OK_REG, true};
  }
};

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

  /*
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

  void reset() {
    _rebuildTimes.assign(_numMaxSamples, 0);
    _currentIndex = 0;
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
  long _rebuildTimesSum = 0;

  /**
   * Stores the current number of samples in the _rebuildTimes vector.
   */
  size_t _numSamples = 0;

  /**
   * Stores the current index in the _rebuildTimes vector where the next sample will be added.
   */
  size_t _currentIndex = 0;

  /**
   * Stores the maximum number of samples in the _rebuildTimes vector.
   * The values is set to 5 to ensure that the estimate is reliable and reflects the current state of the simulation and
   * is not overly influenced by older samples. The magic number 5 was chosen by trial and error with different number
   * of samples.
   */
  size_t _numMaxSamples = 5;
};

/**
 * Estimates a constant value by computing the arithmetic mean of maxN y values with O(N) complexity.
 */
class Mean : public RegressionBase {
 public:
  /**
   * Default constructor is removed as it is mandatory to provide a upper bound on number of sample points.
   */
  Mean() = delete;

  /*
   * Constructor
   */
  Mean(size_t maxN) : RegressionBase(1, maxN) { _y.resize(maxN); }

  void reset() {
    _y.assign(_maxN, 0);
    current_index = 0;
    RegressionBase::reset();
  }

  /**
   *If the maximum number of data points has been reached,no further samples are added to keep the estimator bounded.
   */
  ReturnCode addNewPoint(long const y) {
    _sumY -= _y[current_index];
    _y[current_index] = y;
    _sumY += y;
    current_index = (current_index + 1) % _maxN;
    _n++;

    return ReturnCode::OK_REG;
  }

  /// The prediction of the Mean class takes the mean over the gathered points.
  Result predict() {
    if (not hasEnoughPoints()) {
      return Result{0, ReturnCode::NOT_ENOUGH_POINTS_REG, false};
    }
    auto mean = static_cast<double>(_sumY) / std::min(_maxN, _n);
    return Result{mean, ReturnCode::OK_REG, true};
  }

 private:
  size_t current_index = 0;

  std::vector<long> _y;
};

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
   * Returns the sum of all remainder traversal times, including those with an empty buffer.
   */
  [[nodiscard]] long getSumRemainderTraversalTime() const { return _sumRemainderTraversalTime; }

  /**
   * Reset the manager to the empty state.
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
   */
  [[nodiscard]] double predict(const long numParticlesInBuffer) const {
    if (_numDifferentXConsidered < _minSamples) {
      return 0;
    }
    std::pair<double, double> linRegCoeffs;
    linRegCoeffs = utils::Math::simpleLeastSquares(_numParticlesInBuffer, _remainderTraverlsalTime);
    const double prediction = linRegCoeffs.first + linRegCoeffs.second * static_cast<double>(numParticlesInBuffer);
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
  size_t _currentIndex = 0;

  /**
   * Minimum number of samples needed to make a prediction.
   * Set to 2 as a linear regression requires at least two points to define a line.
   */
  size_t _minSamples = 2;

  /**
   * Maximum number of samples considered to make a prediction.
   * Set to 5 based on trials with different numbers of samples to ensure that the estimate is reliable while keeping
   * computation time low.
   */
  size_t _maxSamples = 5;

  /**
   * Counter to keep track of different number of buffer particles considered for prediction.
   */
  size_t _numDifferentXConsidered = 0;

  /**
   * To keep a running sum of the remainder traversal times for all the iterations, including where the particle buffer
   * size is zero.
   */
  long _sumRemainderTraversalTime = 0;

  /**
   * Used to track the largest _numParticlesInBuffer seen so far.
   * To ensure that incoming _numParticlesInBuffer values are non-decreasing.
   */
  long _currentMaxNumParticlesInBuffer = LONG_MIN;
};

/**
 * Streaming implementation of simple linear regression with x as predictor and y as response variable.
 *
 * Assumes that x values are strictly non-decreasing and positive. This property
 * is exploited to implement a lightweight ring buffer over distinct x values.
 */
class SimpleLinearRegression : public RegressionBase {
 public:
  SimpleLinearRegression() : RegressionBase(2, 5) {
    _y.reserve(_maxN);
    _x.reserve(_maxN);
  }

  /**
   * @param minN Lower limit for the number of points considered for prediction.
   * @param maxN Upper limit for the number of points considered for prediction.
   */
  SimpleLinearRegression(const size_t minN, const size_t maxN) : RegressionBase(minN, maxN) {
    _y.reserve(_maxN);
    _x.reserve(_maxN);
  }

  /// Checks whether at least _minN distinct x values have been gathered, the minimum number required for prediction.
  [[nodiscard]] bool hasEnoughPoints() const { return _numDifferentXConsidered >= _minN; }

  /// Checks whether _maxN distinct x values have been gathered.
  [[nodiscard]] bool reachedMaxPoints() const { return _numDifferentXConsidered >= _maxN; }

  /// Checks whether more than _maxN distinct x values have been gathered.
  [[nodiscard]] bool exceedsMaxPoints() const { return _numDifferentXConsidered > _maxN; }

  /// for tests
  [[nodiscard]] size_t getNumDifferentXConsidered() const { return _numDifferentXConsidered; }

  /// Reset to initial values.
  void reset() {
    RegressionBase::reset();
    _numDifferentXConsidered = 0;
    _currentMaxX = LONG_MIN;
    current_index = 0;
    _x.clear();
    _y.clear();
  }

  /**
   * Maintains a ring buffer over distinct x values.
   *
   * When the buffer is full and a new (larger) x value is added, the oldest
   * (smallest) x/y pair is removed. The accumulated sumY is not treated as a ring buffer.
   */
  ReturnCode addNewPoint(const long x, const long y) {
    // Ensures that added points are monotonically increasing in x, which should be the case anyway.
    if (_currentMaxX < x) {
      _currentMaxX = x;
      /*
       * Ideally, all entries with the same x value should be removed here. However, this case is intentionally
       * ignored, because of performance and code volume, which means there may be more than _maxN values considered
       * for prediction.
       */
      if (not this->reachedMaxPoints()) {
        _x.push_back(x);
        _y.push_back(y);
        _numDifferentXConsidered++;
      } else {
        _x[current_index] = x;
        _y[current_index] = y;
        current_index = (current_index + 1) % _maxN;
      }
    }
    return RegressionBase::addNewPoint(y);
  }

  /**
   * Predicts the y value for a given x using a linear regression over the collected data.
   *
   * @param x Value for which the y value should be returned.
   * @return y A linear regression estimate for x based on the collected data, returned if no error occurred during
   * calculation.
   */
  [[nodiscard]] Result predict(const long x) const {
    if (not this->hasEnoughPoints()) {
      return Result{0, ReturnCode::NOT_ENOUGH_POINTS_REG, false};
    }
    std::pair<double, double> lin_reg_coeffs;
    try {
      lin_reg_coeffs = utils::Math::simpleLeastSquares(_x, _y);
    } catch (...) {
      return Result{0, ReturnCode::UNKNOWN_ERROR_REG, false};
    }
    const double prediction = lin_reg_coeffs.first + lin_reg_coeffs.second * static_cast<double>(x);
    return Result{prediction, ReturnCode::OK_REG, true};
  }

 private:
  /// counter of amount of different x values considered could be actually larger
  size_t _numDifferentXConsidered = 0;

  /**
   * Tracks the largest x seen so far.
   * Assumes that incoming x values are non-decreasing.
   * This allows detecting new distinct x values in O(1).
   */
  long _currentMaxX = LONG_MIN;

  /// Vectors storing the currently considered x values and their corresponding y values for prediction.
  std::vector<long> _x, _y;

  size_t current_index = 0;
};

/**
 * In RebuildDecisionContext, remainder traversals with an empty buffer are ignored for prediction
 * because they could skew the result, but they are still needed for the rebuild decision.
 *
 * This class provides a structure to separate remainder traversal time values with empty buffers from non-zero ones.
 */
class SimpleLinearRegressionSeparateZero {
 public:
  /// Reset both storage structures to their initial values.
  void reset() {
    _zero.reset();
    _rest.reset();
  }

  /// Adds a new point. If x equals zero, store in _zero. Otherwise, store in _rest.
  RegressionBase::ReturnCode addNewPoint(const long x, const long y) {
    RegressionBase::ReturnCode returnCode;
    if (x == 0) {
      returnCode = _zero.addNewPoint(y);
    } else {
      returnCode = _rest.addNewPoint(x, y);
    }
    return returnCode;
  }

  /// For prediction via linear regression, only values from non-empty buffers are considered.
  [[nodiscard]] RegressionBase::Result predict(const long x) { return _rest.predict(x); }

  /**
   * For the overall rebuild decision, the y values from both empty and non-empty buffered remainder traversals are
   * needed.
   */
  [[nodiscard]] long getSumY() const { return _zero.getSumY() + _rest.getSumY(); }

 private:
  /// Stores remainder traversal values with an empty buffer.
  RunningMean _zero;

  /**
   * Stores remainder traversal values with a filled buffer
   * and predicts remainder traversal time for a given buffer filling.
   */
  SimpleLinearRegression _rest{};
};

/**
 * Includes all methods and variables needed to make a rebuild decision by minimizing
 * the sum of rebuild and remainder traversal time spans over the simulation time.
 */
class RebuildDecisionContext {
 public:
  /**
   * Default constructor.
   */
  RebuildDecisionContext() : _rebuildTimeEstimator(), _remainderTraversalTimePredictor(){};

  /// For rebuildNeighborTimeEstimate logging.
  [[nodiscard]] double getRebuildNeighborTimeEstimate() const { return _rebuildNeighborTimeEstimate; }
  /// For remainderTraversalTimeEstimate  logging.
  [[nodiscard]] double getRemainderTraversalTimeEstimate() const { return _remainderTraversalTimeEstimate; }
  /// For particlesBufferEstimate  logging.
  [[nodiscard]] size_t getNumParticlesBufferEstimate() const { return _numParticlesBufferEstimate; }

  /**
   * Called after a neighbor list rebuild to reset remainder traversal estimator state and
   * register the latest rebuild time.
   *
   * @param rebuildTime Measured remainder traversal time
   */
  void afterRebuild(const long rebuildTime) {
    _rebuildTimeEstimator.addNewPoint(rebuildTime);
    _remainderTraversalTimePredictor.reset();
    _remainderTraversalTimeEstimate = std::numeric_limits<double>::quiet_NaN();
  }

  /**
   * The function is called the remainder traversal have been executed for all the functors.
   *
   * Updates the traversal time estimators depending on the current particle
   * buffer size and detects numerical overflow conditions that may trigger
   * a dynamic rebuild.
   *
   * Additionally, tracks the increase in the number of particles in the buffer
   * to support future buffer size estimation.
   *
   * Special case:
   * If a rebuild is performed in the current iteration, the particle buffer
   * is reset due to particles being resorted into containers. In this case,
   * the buffer increase cannot be computed from the immediately preceding
   * iteration and must instead refer to an earlier one.
   *
   * @param remainderTraversalTime Measured total remainder traversal time (over all computeInteraction calls, e.g.
   * 2B+3B)
   * @param numParticlesBuffer     Current number of particles in the particle buffer
   */
  void afterRemainderTraversal(const long remainderTraversalTime, const size_t numParticlesBuffer) {
    _lastRemainderTraversalTime = remainderTraversalTime;

    if (numParticlesBuffer > 0) {
      // numParticleBuffer should not never really decrease, but this is a safety check to avoid underflow in case of
      // unexpected behavior.
      _lastIncreaseNumParticlesBuffer = numParticlesBuffer > _afterFastAddNumParticlesBuffer
                                            ? numParticlesBuffer - _afterFastAddNumParticlesBuffer
                                            : 0;
    }

    _remainderTraversalTimePredictor.addNewPoint(static_cast<long>(numParticlesBuffer), remainderTraversalTime);
  }

  void updateNumParticlesBufferEstimate(const long currentNumParticlesBuffer) {
    _afterFastAddNumParticlesBuffer = currentNumParticlesBuffer;
    _numParticlesBufferEstimate = _afterFastAddNumParticlesBuffer + _lastIncreaseNumParticlesBuffer;
  }

  /**
   * Determines whether the particle buffer has grown large enough such that
   * performing a rebuild is more advantageous than executing another remainder traversal.
   *
   * Updates the internal dynamic rebuild flag (_doDynamicRebuild).
   *
   * @param curr_iter current iteration, or current rebuild frequency if one would rebuild in the current iteration
   * @return True if a dynamic rebuild should be triggered due to particle buffer fullness.
   */
  bool decideToRebuildOnParticleBufferFullness(const unsigned int currIter) {
    bool doDynamicRebuild{false};
    const auto valueRebuild = _rebuildTimeEstimator.predict();

    const auto valueRemainder =
        _remainderTraversalTimePredictor.predict(static_cast<long>(_numParticlesBufferEstimate));

    if (valueRemainder > 0) {
      // for logging
      _rebuildNeighborTimeEstimate = valueRebuild;
      _remainderTraversalTimeEstimate = valueRemainder;

      // First full term: value_rebuild / (curr_iter * (curr_iter + 1)), but denominator is the same, so can be ignored.
      // const double rebuildIncline = value_rebuild / (curr_iter * (curr_iter + 1));

      // Second full term: (value_remainder * curr_iter -
      // static_cast<double>(_remainderTraversalTimePredictor.getSumY())) / (curr_iter * (curr_iter + 1)), but the
      // denominator is same and can be ignored.
      const double remainderIncline =
          (valueRemainder * currIter -
           static_cast<double>(_remainderTraversalTimePredictor.getSumRemainderTraversalTime()));
      if (remainderIncline >= valueRebuild) {
        doDynamicRebuild = true;
      }
    }
    return doDynamicRebuild;
  }

 private:
  /**
   * Handles storing past rebuild times samples and predicting next rebuild time based on simple arithmetic mean of the
   * collected samples.
   */
  RebuildEstimateManager _rebuildTimeEstimator;

  /**
   * Stores (numParticlesBuffer, remainderTraversalTime) pairs for each iteration.
   *
   * Used to predict the remainder traversal time for non-zero buffer sizes,
   * assuming a linear relationship between the particle buffer size (predictor)
   * and the remainder traversal time (response).
   *
   * Has a collector for remainder traversal times for iterations where the particle buffer
   * size is zero. Zero buffer sizes distort linear regression and are therefore
   * handled separately using a mean estimator, because the sum of all times since the last rebuild are still needed.
   */
  RemainderTraversalEstimateManager _remainderTraversalTimePredictor;

  /**
   * Stores the remainder traversal time of the previous iteration.
   * Used by the first decision criterion to estimate the time increase
   * if no rebuild is performed.
   */
  long _lastRemainderTraversalTime{0};

  /// Number of particles in the particle buffer after the fast particles have been added.
  size_t _afterFastAddNumParticlesBuffer{0};

  /**
   * Stores the number of particles added to the buffer during a iteration,
   * due to migrating particles, for estimating the particle buffer filling for the next iteration.
   */
  size_t _lastIncreaseNumParticlesBuffer{0};

  /**
   * Estimated number of particles for the particle buffer based on the last increases
   * and the number of particles in the buffer after fast particles have been added.
   */
  size_t _numParticlesBufferEstimate{0};

  /// For rebuildNeighborTimeEstimate logging.
  double _rebuildNeighborTimeEstimate{std::numeric_limits<double>::quiet_NaN()};
  /// For remainderTraversalTimeEstimate  logging.
  double _remainderTraversalTimeEstimate{std::numeric_limits<double>::quiet_NaN()};
};
}  // namespace autopas::utils::RegressionModels
