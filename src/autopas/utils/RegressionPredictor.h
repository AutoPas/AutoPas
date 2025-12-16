/**
 * @file RegressionPredictor.h
 * @author Natallia Padalinskaya
 * @date 02.12.2025
 */

#pragma once

#include <boost/math/statistics/linear_regression.hpp>

#include "autopas/utils/Math.h"

//#ifdef AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER
namespace autopas::utils {

/**
 * Implements multiple linear regression with two predictor variables.
 * And average of one variable.
 */
class RegressionObject {
 public:
  RegressionObject(const size_t minN, const size_t maxN) : _minN(minN), _maxN(maxN) {}
  enum class ReturnCode {
    OK,
    NOT_ENOUGH_POINTS,
    EXCEEDED_MAX_POINTS,
    OVERFLOW,
    DIVIDE_BY_ZERO,
    UNKNOWN_ERROR,
    REGRESSION_OBJECT
  };
  struct Result {
    double _value;
    ReturnCode _returnCode;
    bool _isOk;
  };
  void reset();
  ReturnCode addNewPoint();
  Result predict();

  [[nodiscard]] bool hasEnoughPoints() const { return _n >= _minN; }
  [[nodiscard]] bool exceedsMaximumPoints() const { return _n >= _maxN; }
  [[nodiscard]] bool isWithinPointsLimits() const { return _n <= _maxN && _n >= _minN; }
  [[nodiscard]] size_t getN() const { return _n; }
  [[nodiscard]] long getSumY() const { return _sumY; }
  void setMinMax(const size_t minN, const size_t maxN) {
    _minN = minN;
    _maxN = maxN;
  }
  void setMaxN(const size_t maxN) { _maxN = maxN; }

 protected:
  size_t _n = 0;
  long _sumY = 0;
  size_t _minN = 0;
  size_t _maxN = 0;
};

class Mean : public RegressionObject {
 public:
  Mean() : RegressionObject(1, 10000) {}

  explicit Mean(const size_t maxN) : RegressionObject(1, maxN) {}

  Mean(const size_t minN, const size_t maxN) : RegressionObject(minN, maxN) {}

  void reset() {
    _sumY = 0;
    _n = 0;
  }

  ReturnCode addNewPoint(long const y) {
    if (exceedsMaximumPoints()) {
      return ReturnCode::EXCEEDED_MAX_POINTS;
    }

    _sumY = Math::safeAdd(_sumY, y);
    if (_sumY == std::numeric_limits<long>::max() || _sumY == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    _n++;
    return ReturnCode::OK;
  }

  Result predict() {
    if (!hasEnoughPoints()) {
      return Result{0, ReturnCode::NOT_ENOUGH_POINTS, false};
    }

    if (!exceedsMaximumPoints()) {
      _lastMean = static_cast<double>(_sumY) / static_cast<double>(_n);
    }
    return Result{_lastMean, ReturnCode::OK, true};
  }

 private:
  double _lastMean = 0.0;
};

class LinearRegression1Predictor : public RegressionObject {
 public:
  LinearRegression1Predictor() : RegressionObject(4, 15000) {}

  LinearRegression1Predictor(const size_t minN, const size_t maxN) : RegressionObject(minN, maxN) {}

  void reset() {
    _sumX = 0.0;
    _sumY = 0.0;
    _sumXX = 0.0;
    _sumXY = 0.0;
    _n = 0;
  }
  ReturnCode addNewPoint(const long x, const long y) {
    if (exceedsMaximumPoints()) {
      return ReturnCode::EXCEEDED_MAX_POINTS;
    }

    _sumX = Math::safeAdd(_sumX, x);
    if (_sumX == std::numeric_limits<long>::max() || _sumX == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    _sumY = Math::safeAdd(_sumY, y);
    if (_sumY == std::numeric_limits<long>::max() || _sumY == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    long const x_squared = Math::safeMul(x, x);
    if (x_squared == std::numeric_limits<long>::max() || x_squared == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    _sumXX = Math::safeAdd(_sumXX, x_squared);
    if (_sumXX == std::numeric_limits<long>::max() || _sumXX == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    long xy = Math::safeMul(x, y);
    if (xy == std::numeric_limits<long>::max() || xy == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    _sumXY = Math::safeAdd(_sumXY, xy);
    if (_sumXY == std::numeric_limits<long>::max() || _sumXY == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    _n++;
    return ReturnCode::OK;
  }

  Result predict(const long x) {
    if (!hasEnoughPoints()) {
      return Result{0, ReturnCode::NOT_ENOUGH_POINTS, false};
    }
    if (!exceedsMaximumPoints()) {
      // Calculate denominator: n * sumXX - sumX * sumX
      const long n_long = static_cast<long>(_n);

      const long n_sumXX = Math::safeMul(n_long, _sumXX);
      if (n_sumXX == std::numeric_limits<long>::max() || n_sumXX == std::numeric_limits<long>::min()) {
        return Result{0, ReturnCode::OVERFLOW, false};
      }

      const long sumX_squared = Math::safeMul(_sumX, _sumX);
      if (sumX_squared == std::numeric_limits<long>::max() || sumX_squared == std::numeric_limits<long>::min()) {
        return Result{0, ReturnCode::OVERFLOW, false};
      }

      const long denominator = Math::safeSub(n_sumXX, sumX_squared);
      if (denominator == std::numeric_limits<long>::max() || denominator == std::numeric_limits<long>::min()) {
        return Result{0, ReturnCode::OVERFLOW, false};
      }

      if (denominator == 0) {
        return Result{0, ReturnCode::DIVIDE_BY_ZERO, false};
      }

      // Calculate numerator: n * sumXY - sumX * sumY
      const long n_sumXY = Math::safeMul(n_long, _sumXY);
      if (n_sumXY == std::numeric_limits<long>::max() || n_sumXY == std::numeric_limits<long>::min()) {
        return Result{0, ReturnCode::OVERFLOW, false};
      }

      const long sumX_sumY = Math::safeMul(_sumX, _sumY);
      if (sumX_sumY == std::numeric_limits<long>::max() || sumX_sumY == std::numeric_limits<long>::min()) {
        return Result{0, ReturnCode::OVERFLOW, false};
      }

      const long numerator = Math::safeSub(n_sumXY, sumX_sumY);
      if (numerator == std::numeric_limits<long>::max() || numerator == std::numeric_limits<long>::min()) {
        return Result{0, ReturnCode::OVERFLOW, false};
      }

      // Calculate beta1 and beta0
      _last_ß1 = static_cast<double>(numerator) / static_cast<double>(denominator);
      _last_ß0 = static_cast<double>(_sumY) / static_cast<double>(_n) -
                 _last_ß1 * (static_cast<double>(_sumX) / static_cast<double>(_n));
    }
    const double prediction = _last_ß0 + _last_ß1 * static_cast<double>(x);
    return Result{prediction, ReturnCode::OK, true};
  }

 private:
  long _sumX = 0;
  long _sumXX = 0;
  long _sumXY = 0;
  double _last_ß0 = 0.0;
  double _last_ß1 = 0.0;
};

class LinearRegression1PredictorBoost : public RegressionObject {
 public:
  LinearRegression1PredictorBoost() : RegressionObject(2, 100) {
    _y.reserve(_maxN);
    _x.reserve(_maxN);
  }

  LinearRegression1PredictorBoost(const size_t minN, const size_t maxN) : RegressionObject(minN, maxN) {
    _y.reserve(_maxN);
    _x.reserve(_maxN);
  }

  [[nodiscard]] bool hasEnoughPoints() const { return _numDifferentX >= _minN; }
  [[nodiscard]] bool exceedsMaximumPoints() const { return _numDifferentX >= _maxN; }
  [[nodiscard]] bool isWithinPointsLimits() const { return _numDifferentX <= _maxN && _numDifferentX >= _minN; }

  void reset() {
    _numDifferentX = 0;
    _currentMaxX = LONG_MIN;
    _sumY = 0;
    _n = 0;
    _x.clear();
    _y.clear();
  }
  ReturnCode addNewPoint(const long x, const long y) {
    _sumY = Math::safeAdd(_sumY, y);
    if (_sumY == std::numeric_limits<long>::max() || _sumY == std::numeric_limits<long>::min()) {
      return ReturnCode::OVERFLOW;
    }

    if (_currentMaxX < x) {
      _numDifferentX++;
      _currentMaxX = x;
      if (this->exceedsMaximumPoints()) {
        _x.erase(_x.begin());
        _y.erase(_y.begin());
      }
    }

    _x.push_back(x);
    _y.push_back(y);

    _n++;
    return ReturnCode::OK;
  }

  Result predict(const long x) {
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
  size_t _numDifferentX = 0;
  long _currentMaxX = LONG_MIN;
  std::vector<long> _x;
  std::vector<long> _y;
};

}  // namespace autopas::utils
//#endif
