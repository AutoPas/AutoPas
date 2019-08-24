/**
 * @file AcquisitionFunctionOption.h
 * @author Jan Nguyen
 * @date 26.06.19
 */

#pragma once

#include <set>

namespace autopas {

/**
 * Different acquisition functions
 */
enum AcquisitionFunctionOption {
  upperConfidenceBound,
  lowerConfidenceBound,
  mean,
  variance,
  probabilityOfDecrease,
  expectedDecrease
};

/**
 * Provides a way to iterate over the possible choices of AcquisitionFunction.
 */
static const std::set<AcquisitionFunctionOption> allAcquisitionFunctionOptions = {
    AcquisitionFunctionOption::upperConfidenceBound,
    AcquisitionFunctionOption::lowerConfidenceBound,
    AcquisitionFunctionOption::mean,
    AcquisitionFunctionOption::variance,
    AcquisitionFunctionOption::probabilityOfDecrease,
    AcquisitionFunctionOption::expectedDecrease};
}  // namespace autopas
