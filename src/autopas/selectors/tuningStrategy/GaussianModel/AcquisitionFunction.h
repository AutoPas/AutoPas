/**
 * @file AcquisitionFunction.h
 * @author Jan Nguyen
 * @date 13.04.20
 */

#pragma once

#include "autopas/options/AcquisitionFunctionOption.h"

/**
 * Functions used for acquisition prediction used in Gaussian Process models
 */
namespace autopas::AcquisitionFunction {
/**
 * Calculate given acquisition function
 * @param af acquisition function
 * @param mean
 * @param var
 * @param currentMax highest output-value in the current evidence set
 * @return
 */
double calcAcquisition(const AcquisitionFunctionOption &af, double mean, double var, double currentMax = 0.);
}  // namespace autopas::AcquisitionFunction
