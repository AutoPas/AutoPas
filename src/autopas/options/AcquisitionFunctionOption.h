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
  UpperConfidenceBound,
  LowerConfidenceBound,
  Mean,
  Variance,
  ProbabilityOfDecrease,
  ExpectedDecrease
};

/**
 * Provides a way to iterate over the possible choices of AcquisitionFunction.
 */
static const std::set<AcquisitionFunctionOption> allAcquisitionFunctionOptions = {
    AcquisitionFunctionOption::UpperConfidenceBound,
    AcquisitionFunctionOption::LowerConfidenceBound,
    AcquisitionFunctionOption::Mean,
    AcquisitionFunctionOption::Variance,
    AcquisitionFunctionOption::ProbabilityOfDecrease,
    AcquisitionFunctionOption::ExpectedDecrease};
}  // namespace autopas
