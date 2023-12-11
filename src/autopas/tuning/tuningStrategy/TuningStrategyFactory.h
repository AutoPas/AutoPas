/**
 * @file TuningStrategyFactory.h
 * @author seckler
 * @date 07.02.2020
 */

#pragma once

#include <memory>
#include <set>
#include <string>

#include "autopas/options/TuningStrategyOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyFactoryInfo.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas::TuningStrategyFactory {
/**
 * Generates a new Tuning Strategy object.
 * @param searchSpace
 * @param tuningStrategyOption
 * @param info
 * @param outputSuffix
 * @return Pointer to the tuning strategy object or the null pointer if an exception was suppressed.
 */
std::unique_ptr<TuningStrategyInterface> generateTuningStrategy(const std::set<Configuration> &searchSpace,
                                                                TuningStrategyOption tuningStrategyOption,
                                                                const TuningStrategyFactoryInfo &info,
                                                                const std::string &outputSuffix = "");
}  // namespace autopas::TuningStrategyFactory
