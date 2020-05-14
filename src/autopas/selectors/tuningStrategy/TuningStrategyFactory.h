/**
 * @file TuningStrategyFactory.h
 * @author seckler
 * @date 07.02.20
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/utils/WrapMPI.h"

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
    autopas::TuningStrategyOption tuningStrategyOption, const std::set<autopas::ContainerOption> &allowedContainers,
    autopas::NumberSet<double> &allowedCellSizeFactors, const std::set<autopas::TraversalOption> &allowedTraversals,
    const std::set<autopas::DataLayoutOption> &allowedDataLayouts,
    const std::set<autopas::Newton3Option> &allowedNewton3Options, unsigned int maxEvidence,
    AcquisitionFunctionOption acquisitionFunctionOption, MPI_Comm);
}  // namespace autopas::TuningStrategyFactory
