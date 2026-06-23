/**
 * @file TuningStrategyFactory.h
 * @author seckler
 * @date 07.02.2020
 */

#pragma once

#include "TuningStrategyFactoryInfo.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"

namespace autopas::TuningStrategyFactory {
/**
 * Generates a new Tuning Strategy object.
 * @param searchSpace Search space of algorithmic configurations to choose from.
 * @param tuningStrategyOption the tuning strategy type.
 * @param info TuningStrategyFactoryInfo containing information which may be relevant when construction the given
 * tuning strategy.
 * @param interactionType Type of interaction.
 * @param outputSuffix
 * @return Pointer to the tuning strategy object or the null pointer if an exception was suppressed.
 */
std::unique_ptr<TuningStrategyInterface> generateTuningStrategy(const std::set<Configuration> &searchSpace,
                                                                TuningStrategyOption tuningStrategyOption,
                                                                const TuningStrategyFactoryInfo &info,
                                                                InteractionTypeOption interactionType,
                                                                const std::string &outputSuffix = "");
}  // namespace autopas::TuningStrategyFactory
