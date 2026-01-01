/**
 * @file RegressionModels.h
 * @author Natallia Padalinskaya
 * @date 02.12.2025
 */

#pragma once

#include <boost/math/statistics/linear_regression.hpp>

#include "autopas/utils/Math.h"

namespace autopas::utils {

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
  // TODO: minN < maxN
  RegressionBase(const size_t minN, const size_t maxN) : _minN(minN), _maxN(maxN) {}
  enum class ReturnCode { OK, NOT_ENOUGH_POINTS, EXCEEDED_MAX_POINTS, OVERFLOW, DIVIDE_BY_ZERO, UNKNOWN_ERROR };
  struct Result {
    double _value;
    ReturnCode _returnCode;
    bool _isOk;
  };
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
      return ReturnCode::OVERFLOW;
    }
    incrementN();
    return ReturnCode::OK;
  }

  void incrementN() { _n++; }
  static bool isOk(const ReturnCode return_code) {
    if (return_code == ReturnCode::OK) {
      return true;
    }
    return false;
  }

  /// Checks whether a value indicates an arithmetic overflow based on sentinel limits.
  static bool isOverflown(const long x) {
    if (x == std::numeric_limits<long>::max() || x == std::numeric_limits<long>::min()) {
      return true;
    }
    return false;
  }

  [[nodiscard]] bool hasEnoughPoints() const { return _n >= _minN; }
  [[nodiscard]] bool reachedMaxPoints() const { return _n >= _maxN; }
  [[nodiscard]] bool exceedsMaxPoints() const { return _n > _maxN; }
  [[nodiscard]] bool isWithinPointsLimits() const { return _n <= _maxN && _n >= _minN; }
  [[nodiscard]] size_t getN() const { return _n; }
  [[nodiscard]] long getSumY() const { return _sumY; }

  // for tuning minN and maxN through YAML parameters
  void setMinMax(const size_t minN, const size_t maxN) {
    _minN = minN;
    _maxN = maxN;
  }
  void setMaxN(const size_t maxN) { _maxN = maxN; }

 protected:
  /// Number of data points currently stored
  size_t _n = 0;

  long _sumY = 0;

  /// Amount of considered points for prediction has a lower (_minN) and an upper (_maxN) limit.
  size_t _minN = 0;
  size_t _maxN = 0;
};

/**
 * Estimates a constant value by computing the arithmetic mean of observed
 * y values in a streaming fashion.
 *
 * Used as a lightweight baseline estimator for quantities with low variance
 * (e.g. rebuildNeighborListTime).
 */
class Mean : public RegressionBase {
 public:
  // TODO: minN > 0
  Mean() : RegressionBase(1, 100) {}

  explicit Mean(const size_t maxN) : RegressionBase(1, maxN) {}

  Mean(const size_t minN, const size_t maxN) : RegressionBase(minN, maxN) {}

  void reset() { RegressionBase::reset(); }

  /**
   *If the maximum number of data points has been reached,no further samples are added to keep the estimator bounded.
   */
  ReturnCode addNewPoint(long const y) {
    if (reachedMaxPoints()) {
      return ReturnCode::EXCEEDED_MAX_POINTS;
    }

    return RegressionBase::addNewPoint(y);
  }

  Result predict() {
    if (!hasEnoughPoints()) {
      return Result{0, ReturnCode::NOT_ENOUGH_POINTS, false};
    }

    // If the maximum number of points has not been reached yet,
    // recompute the mean using the current accumulated values.
    // Once the estimator is full, the cached mean can be reused.
    if (!exceedsMaxPoints()) {
      _lastMean = static_cast<double>(_sumY) / static_cast<double>(_n);
    }
    return Result{_lastMean, ReturnCode::OK, true};
  }

 private:
  /// Cached mean value to avoid recomputation once the maximum number of points is reached
  double _lastMean = 0.0;
};

/**
 * Streaming implementation of simple linear regression using the Boost library with x as predictor
 * and y as response variable.
 *
 * Assumes that x values are strictly non-decreasing and positive. This property
 * is exploited to implement a lightweight ring buffer over distinct x values.
 */
class SimpleLinearRegressionBoost : public RegressionBase {
 public:
  SimpleLinearRegressionBoost() : RegressionBase(2, 5) {
    _y.reserve(_maxN);
    _x.reserve(_maxN);
  }

  SimpleLinearRegressionBoost(const size_t minN, const size_t maxN) : RegressionBase(minN, maxN) {
    _y.reserve(_maxN);
    _x.reserve(_maxN);
  }

  /// Point limits are evaluated based on the number of distinct x values,
  /// not on the total number of samples.
  [[nodiscard]] bool hasEnoughPoints() const { return _numDifferentXConsidered >= _minN; }
  [[nodiscard]] bool reachedMaxPoints() const { return _numDifferentXConsidered >= _maxN; }
  [[nodiscard]] bool exceedsMaxPoints() const { return _numDifferentXConsidered > _maxN; }
  [[nodiscard]] bool isWithinPointsLimits() const {
    return _numDifferentXConsidered <= _maxN && _numDifferentXConsidered >= _minN;
  }

  // for tests
  [[nodiscard]] size_t getNumDifferentXConsidered() const { return _numDifferentXConsidered; }

  void reset() {
    RegressionBase::reset();
    _numDifferentXConsidered = 0;
    _currentMaxX = LONG_MIN;
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
    if (_currentMaxX < x) {
      _currentMaxX = x;
      if (this->reachedMaxPoints()) {
        // Ideally, all entries with the same x value should be removed here. However, this case is intentionally
        // ignored, because of performance and code volume, which means there may be more than _maxN values considered
        // for prediction.
        _x.erase(_x.begin());
        _y.erase(_y.begin());
      } else {
        _numDifferentXConsidered++;
      }
    }

    _x.push_back(x);
    _y.push_back(y);

    return RegressionBase::addNewPoint(y);
  }

  [[nodiscard]] Result predict(const long x) const {
    if (not this->hasEnoughPoints()) {
      return Result{0, ReturnCode::NOT_ENOUGH_POINTS, false};
    }
    std::pair<double, double> ßs;
    try {
      ßs = boost::math::statistics::simple_ordinary_least_squares(_x, _y);
    } catch (...) {
      return Result{0, ReturnCode::UNKNOWN_ERROR, false};
    }
    const double prediction = ßs.first + ßs.second * static_cast<double>(x);
    return Result{prediction, ReturnCode::OK, true};
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

  std::vector<long> _x;
  std::vector<long> _y;
};

class SimpleLinearRegressionSeparateZero {
 public:
  void reset() {
    _zero.reset();
    _rest.reset();
  }

  RegressionBase::ReturnCode addNewPoint(const long x, const long y) {
    RegressionBase::ReturnCode returnCode;
    if (x == 0) {
      returnCode = _zero.addNewPoint(y);
    } else {
      returnCode = _rest.addNewPoint(x, y);
    }
    return returnCode;
  }

  [[nodiscard]] RegressionBase::Result predict(const long x) { return _rest.predict(x); }

  [[nodiscard]] long getSumY() const { return _zero.getSumY() + _rest.getSumY(); }

  void setMinMax(const size_t minN, const size_t maxN) { _rest.setMinMax(minN, maxN); }

 private:
  Mean _zero{SIZE_MAX};
  SimpleLinearRegressionBoost _rest{};
};

class RebuildDecisionContext {
 public:
  // log
  [[nodiscard]] double getRebuildNeighborTimeEstimate() const { return _rebuildNeighborTimeEstimate; }
  [[nodiscard]] double getRemainderTraversalTimeEstimate() const { return _remainderTraversalTimeEstimate; }
  [[nodiscard]] size_t getNumParticlesBufferEstimate() const { return _numParticlesBufferEstimate; }

  /**
   * Called after a neighbor list rebuild to reset remainder traversal estimator state and
   * register the latest rebuild time.
   *
   * @param rebuildTime Measured remainder traversal time
   * @param doTuningRebuild
   * @param isFirstIteration
   */
  void afterRebuild(const long rebuildTime, const bool doTuningRebuild, const bool isFirstIteration) {
    // tuning iterations distort the rebuild time estimate
    if (!doTuningRebuild and !isFirstIteration) {
      if (_rebuildNeighborTimeMean.addNewPoint(rebuildTime) == RegressionBase::ReturnCode::OVERFLOW) {
        _rebuildNeighborTimeMean.reset();
        _rebuildNeighborTimeEstimate = std::numeric_limits<double>::quiet_NaN();
      }
    }

    _remainderTraversalTimePredictor.reset();

    _remainderTraversalTimeEstimate = std::numeric_limits<double>::quiet_NaN();
  }

  /**
   * Called after the remainder traversal has been executed.
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
   * @param remainderTraversalTime Measured remainder traversal time
   * @param numParticlesBuffer     Current number of particles in the particle buffer
   *
   * @return True if a dynamic rebuild should be triggered in the next iteration (e.g. due to overflow)
   */
  bool afterRemainderTraversal(const long remainderTraversalTime, const size_t numParticlesBuffer) {
    _lastRemainderTraversalTime = remainderTraversalTime;

    if (numParticlesBuffer > 0) {
      _lastIncreaseNumParticlesBuffer = numParticlesBuffer - _afterFastAddNumParticlesBuffer;
    }

    return _remainderTraversalTimePredictor.addNewPoint(static_cast<long>(numParticlesBuffer),
                                                        remainderTraversalTime) == RegressionBase::ReturnCode::OVERFLOW;
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
   * @param rf rebuild frequency, if one would rebuild in the current iteration
   * @return True if a dynamic rebuild should be triggered due to particle buffer fullness,if an initial rebuild time
   * estimate is required, or if an error or overflow occurred.
   */
  bool decideToRebuildOnParticleBufferFullness(const unsigned int rf) {
    bool doDynamicRebuild{false};
    const auto [value_rebuild, returnCode_rebuild, isOk_rebuild] = _rebuildNeighborTimeMean.predict();

    const auto [value_remainder, returnCode_remainder, isOk_remainder] =
        _remainderTraversalTimePredictor.predict(static_cast<long>(_numParticlesBufferEstimate));

    if (isOk_rebuild and isOk_remainder and value_remainder > 0) {
      // log
      _rebuildNeighborTimeEstimate = value_rebuild;
      _remainderTraversalTimeEstimate = value_remainder;

      const double rebuildIncline = value_rebuild / (rf * rf);
      const double remainderIncline =
          (value_remainder - static_cast<double>(_remainderTraversalTimePredictor.getSumY()) / rf) / (rf + 1);
      if (remainderIncline >= rebuildIncline or value_rebuild <= value_remainder) {
        doDynamicRebuild = true;
      }

    } else if (returnCode_rebuild == RegressionBase::ReturnCode::OVERFLOW ||
               returnCode_remainder == RegressionBase::ReturnCode::OVERFLOW ||
               // No rebuild has been performed yet; trigger one to obtain an initial rebuild time estimate
               returnCode_rebuild == RegressionBase::ReturnCode::NOT_ENOUGH_POINTS ||
               // if the boost library linear regression implementation throws an error
               returnCode_remainder == RegressionBase::ReturnCode::UNKNOWN_ERROR) {
      doDynamicRebuild = true;
    }
    return doDynamicRebuild;
  }

 private:
  /**
   * Stores rebuild neighbor list times.
   * The rebuild time is assumed to be approximately constant and is therefore
   * estimated using a running mean.
   */
  Mean _rebuildNeighborTimeMean{};

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
  SimpleLinearRegressionSeparateZero _remainderTraversalTimePredictor{};

  /**
   * Stores the remainder traversal time of the previous iteration.
   * Used by the first decision criterion to estimate the time increase
   * if no rebuild is performed.
   */
  long _lastRemainderTraversalTime{0};

  /**
   * Estimates the particle buffer size by extrapolating from the particle increase
   * between iterations.
   *
   * The rebuild decision must be made before particles migrate between containers,
   * including across periodic boundaries.
   */
  size_t _afterFastAddNumParticlesBuffer{0};
  size_t _lastIncreaseNumParticlesBuffer{0};
  size_t _numParticlesBufferEstimate{0};

  // log
  double _rebuildNeighborTimeEstimate{std::numeric_limits<double>::quiet_NaN()};
  double _remainderTraversalTimeEstimate{std::numeric_limits<double>::quiet_NaN()};
};
}  // namespace autopas::utils
