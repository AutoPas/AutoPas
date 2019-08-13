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
  /**
   * Upper confidence bound
   */
  ucb,
  /**
   * Lower confidence bound
   */
  lcb,
  /**
   * mean
   */
  mean,
  /**
   * variance
   */
  var
};

/**
 * Provides a way to iterate over the possible choices of AcquisitionFunction.
 */
static const std::set<AcquisitionFunctionOption> allAcquisitionFunctionOptions = {
    AcquisitionFunctionOption::ucb, AcquisitionFunctionOption::lcb, AcquisitionFunctionOption::mean,
    AcquisitionFunctionOption::var};
}  // namespace autopas
