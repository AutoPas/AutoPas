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
namespace AcquisitionFunction {
/**
 * Calculate given acquisition function
 * @param af acquisition function
 * @param mean
 * @param var
 * @param currentMax highest output-value in the current evidence set
 * @return
 */
static double calcAcquisition(AcquisitionFunctionOption af, double mean, double var, double currentMax = 0.) {
  switch (af) {
    case AcquisitionFunctionOption::upperConfidenceBound: {
      return mean + 2 * std::sqrt(var);
    }
    case AcquisitionFunctionOption::mean: {
      return mean;
    }
    case AcquisitionFunctionOption::variance: {
      return var;
    }
    case AcquisitionFunctionOption::probabilityOfImprovement: {
      return utils::Math::normalCDF((mean - currentMax) / std::sqrt(var));
    }
    case AcquisitionFunctionOption::expectedImprovement: {
      double stddev = std::sqrt(var);
      double maxNormed = (mean - currentMax) / stddev;
      double cdf = utils::Math::normalCDF(maxNormed);
      double pdf = utils::Math::normalPDF(maxNormed);
      return (mean - currentMax) * cdf + stddev * pdf;
    }
  }

  autopas::utils::ExceptionHandler::exception("AcquisitionFunction.calcAcquisition: Unknown acquisition function {}.",
                                              af);
  return 0;
}
}  // namespace AcquisitionFunction
}  // namespace autopas
