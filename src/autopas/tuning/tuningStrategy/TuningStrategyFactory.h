/**
 * @file TuningStrategyFactory.h
 * @author seckler
 * @date 07.02.2020
 */

#pragma once

#include "TuningStrategyFactoryInfo.h"
#include "TuningStrategyInterface.h"
#include "autopas/options/TuningStrategyOption.h"

namespace autopas::TuningStrategyFactory {
/**
 * Generates a new Tuning Strategy object.
 * @param searchSpace
 * @param tuningStrategyOption
 * @param info
 * @param outputSuffix
 * @return Pointer to the tuning strategy object or the null pointer if an exception was suppressed.
 */
std::unique_ptr<autopas::TuningStrategyInterface> generateTuningStrategy(
    const std::set<Configuration> &searchSpace, autopas::TuningStrategyOption tuningStrategyOption,
    const TuningStrategyFactoryInfo &info, const std::string &outputSuffix = "");
}  // namespace autopas::TuningStrategyFactory
