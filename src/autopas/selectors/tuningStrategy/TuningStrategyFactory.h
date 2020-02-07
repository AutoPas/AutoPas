/**
 * @file TuningStrategyFactory.h
 * @author seckler
 * @date 07.02.20
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/TuningStrategyOption.h"

namespace autopas::TuningStrategyFactory {
/**
 * Generates a new Tuning Strategy object.
 * @param tuningStrategyOption
 * @param allowedContainers
 * @param allowedCellSizeFactors
 * @param allowedTraversals
 * @param allowedDataLayouts
 * @param allowedNewton3Options
 * @param maxEvidence
 * @param acquisitionFunctionOption
 * @return Pointer to the tuning strategy object or the nullpointer if an exception was suppressed.
 */
std::unique_ptr<autopas::TuningStrategyInterface> generateTuningStrategy(
    TuningStrategyOption tuningStrategyOption, const std::set<ContainerOption> &allowedContainers,
    NumberSet<double> &allowedCellSizeFactors, const std::set<TraversalOption> &allowedTraversals,
    const std::set<DataLayoutOption> &allowedDataLayouts, const std::set<Newton3Option> &allowedNewton3Options,
    unsigned int maxEvidence, AcquisitionFunctionOption acquisitionFunctionOption);
}  // namespace autopas::TuningStrategyFactory
