/**
 * @file AcquisitionFunction.h
 * @author Jan Nguyen
 * @date 13.04.20
 */

#pragma once

#include <Eigen/Core>
#include <utility>

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"

namespace autopas {

/**
 * Functions used for acquisition prediction used in Gaussian Process models
 */
class AcquisitionFunction {
 public:
  /**
   * @brief Calculate given acquisition function
   * @param af acquisition function
   * @param mean
   * @param var
   * @param target
   * @return
   */
  static double calcAcquisition(AcquisitionFunctionOption af, double mean, double var, double target = 0.) {
    switch (af) {
      case AcquisitionFunctionOption::upperConfidenceBound: {
        return mean + 2 * std::sqrt(var);
      }
      case AcquisitionFunctionOption::lowerConfidenceBound: {
        return mean - 2 * std::sqrt(var);
      }
      case AcquisitionFunctionOption::mean: {
        return mean;
      }
      case AcquisitionFunctionOption::variance: {
        return var;
      }
      case AcquisitionFunctionOption::probabilityOfDecrease: {
        return utils::Math::normalCDF((target - mean) / std::sqrt(var));
      }
      case AcquisitionFunctionOption::expectedDecrease: {
        double stddev = std::sqrt(var);
        double minNormed = (target - mean) / stddev;
        return (target - mean) * utils::Math::normalCDF(minNormed) + stddev * utils::Math::normalPDF(minNormed);
      }
    }

    autopas::utils::ExceptionHandler::exception("AcquisitionFunction.calcAcquisition: Unknown acquisition function {}.",
                                                af);
    return 0;
  }
};
}  // namespace autopas
