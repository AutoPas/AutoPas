/**
 * @file configureAutoPas.cpp
 * @date 03.11.2026
 */
#include "configureAutoPas.h"

// The main AutoPas class is instantiated in templateInstantiations/AutoPasClass.cpp.
//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
//! @endcond

void configureAutoPas(autopas::AutoPas<ParticleType> &autoPas, const MDFlexConfig &configuration,
                      const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                      const std::string &outputSuffix) {
  autoPas.setAllowedCellSizeFactors(*configuration.cellSizeFactors.value);
  autoPas.setAllowedContainers(configuration.containerOptions.value);
  autoPas.setAllowedInteractionTypeOptions(configuration.getInteractionTypes());
  // Pairwise specific options
  autoPas.setAllowedDataLayouts(configuration.dataLayoutOptions.value,
                                           autopas::InteractionTypeOption::pairwise);
  autoPas.setAllowedNewton3Options(configuration.newton3Options.value,
                                              autopas::InteractionTypeOption::pairwise);
  autoPas.setAllowedTraversals(configuration.traversalOptions.value,
                                          autopas::InteractionTypeOption::pairwise);
  autoPas.setAllowedVecPatterns(configuration.vecPatternOptions.value,
                                           autopas::InteractionTypeOption::pairwise);
  autoPas.setAllowedLoadEstimators(configuration.loadEstimatorOptions.value);
  // Triwise specific options
  autoPas.setAllowedDataLayouts(configuration.dataLayoutOptions3B.value,
                                           autopas::InteractionTypeOption::triwise);
  autoPas.setAllowedNewton3Options(configuration.newton3Options3B.value,
                                              autopas::InteractionTypeOption::triwise);
  autoPas.setAllowedTraversals(configuration.traversalOptions3B.value,
                                          autopas::InteractionTypeOption::triwise);
  // General options
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(configuration.cutoff.value);
  autoPas.setRelativeOptimumRange(configuration.relativeOptimumRange.value);
  autoPas.setMaxTuningPhasesWithoutTest(configuration.maxTuningPhasesWithoutTest.value);
  autoPas.setRelativeBlacklistRange(configuration.relativeBlacklistRange.value);
  autoPas.setEvidenceFirstPrediction(configuration.evidenceFirstPrediction.value);
  autoPas.setExtrapolationMethodOption(configuration.extrapolationMethodOption.value);
  autoPas.setNumSamples(configuration.tuningSamples.value);
  autoPas.setEarlyStoppingFactor(configuration.earlyStoppingFactor.value);
  autoPas.setMaxEvidence(configuration.tuningMaxEvidence.value);
  autoPas.setRuleFileName(configuration.ruleFilename.value);
  autoPas.setFuzzyRuleFileName(configuration.fuzzyRuleFilename.value);
  autoPas.setSelectorStrategy(configuration.selectorStrategy.value);
  autoPas.setTuningInterval(configuration.tuningInterval.value);
  autoPas.setTuningStrategyOption(configuration.tuningStrategyOptions.value);
  autoPas.setTuningMetricOption(configuration.tuningMetricOption.value);
  autoPas.setUseLOESSSmoothening(configuration.useLOESSSmoothening.value);
  autoPas.setEnergySensorOption(configuration.energySensorOption.value);
  autoPas.setMPITuningMaxDifferenceForBucket(configuration.MPITuningMaxDifferenceForBucket.value);
  autoPas.setMPITuningWeightForMaxDensity(configuration.MPITuningWeightForMaxDensity.value);
  autoPas.setVerletClusterSize(configuration.verletClusterSize.value);
  autoPas.setVerletRebuildFrequency(configuration.verletRebuildFrequency.value);
  autoPas.setVerletSkin(configuration.verletSkinRadius.value);
  autoPas.setDeltaT(configuration.deltaT.value);
  autoPas.setUseTuningLogger(configuration.useTuningLogger.value);
  autoPas.setAoSSortingThreshold(configuration.aosSortingThreshold.value);
  autoPas.setSoASortingThreshold(configuration.soaSortingThreshold.value);
  autoPas.setUseSortingThresholdBenchmark(configuration.useSortingThresholdBenchmark.value);
  autoPas.setOutputSuffix(outputSuffix);
  autopas::Logger::get()->set_level(configuration.logLevel.value);
}
