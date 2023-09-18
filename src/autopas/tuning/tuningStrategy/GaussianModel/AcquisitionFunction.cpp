/**
 * @file AcquisitionFunction.cpp
 * @author F. Gratl
 * @date 17.11.22
 */

#include "AcquisitionFunction.h"

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"

double autopas::AcquisitionFunction::calcAcquisition(const AcquisitionFunctionOption &af, double mean, double var,
                                                     double currentMax) {
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
      const double stddev = std::sqrt(var);
      const double maxNormed = (mean - currentMax) / stddev;
      const double cdf = utils::Math::normalCDF(maxNormed);
      const double pdf = utils::Math::normalPDF(maxNormed);
      return (mean - currentMax) * cdf + stddev * pdf;
    }
  }

  autopas::utils::ExceptionHandler::exception("AcquisitionFunction.calcAcquisition: Unknown acquisition function {}.",
                                              af);
  return 0;
}
