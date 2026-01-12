/**
 * @file BayesianSearch.cpp
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "BayesianSearch.h"

#include <algorithm>
#include <iterator>

#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/utils/StringUtils.h"

autopas::BayesianSearch::BayesianSearch(
    const InteractionTypeOption &interactionType, const std::set<ContainerOption> &allowedContainerOptions,
    const autopas::NumberSet<double> &allowedCellSizeFactors, const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options,
    const std::set<VectorizationPatternOption> &allowedVecPatternOptions, size_t maxEvidence,
    autopas::AcquisitionFunctionOption predAcqFunction, size_t predNumLHSamples, unsigned long seed)
    : _interactionType(interactionType),
      _containerOptionsSet(allowedContainerOptions),
      _dataLayoutOptions(allowedDataLayoutOptions.begin(), allowedDataLayoutOptions.end()),
      _vecPatternOptions(allowedVecPatternOptions.begin(), allowedVecPatternOptions.end()),
      _newton3Options(allowedNewton3Options.begin(), allowedNewton3Options.end()),
      _cellSizeFactors(allowedCellSizeFactors.clone()),
      _encoder(),
      _invalidConfigs(),
      _rng(seed),
      _gaussianProcess(0, 0.01, _rng),
      _maxEvidence(maxEvidence),
      _predAcqFunction(predAcqFunction),
      _predNumLHSamples(predNumLHSamples) {
  if (predNumLHSamples <= 0) {
    utils::ExceptionHandler::exception(
        "BayesianSearch: Number of samples used for predictions must be greater than 0!");
  }

  for (const auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption, _interactionType);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (const auto &traversalOption : allowedAndApplicable) {
      // if load estimators are not applicable LoadEstimatorOption::none is returned.
      const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
          loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
      for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
        _containerTraversalEstimatorOptions.emplace_back(containerOption, traversalOption, loadEstimatorOption);
      }
    }
  }

  if (searchSpaceIsEmpty()) {
    autopas::utils::ExceptionHandler::exception("BayesianSearch: No valid configurations could be created.");
  }

  _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options,
                             *_cellSizeFactors, _vecPatternOptions);
  _gaussianProcess.setDimension(_encoder.getOneHotDims());
}

bool autopas::BayesianSearch::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                  const EvidenceCollection &evidenceCollection) {
  // if enough evidence was collected abort the tuning process.
  if (_gaussianProcess.numEvidence() >= _maxEvidence) {
    configQueue.clear();
    return true;
  }

  // Sample the search space, check that the samples are in the available configurations, and
  // sort the remaining intersection of the config queue
  for (size_t i = 0; i < maxAttempts; ++i) {
    auto currentSamples = sampleAcquisitions(_predNumLHSamples, _predAcqFunction);

    // Filter out all configurations which are either:
    //   - invalid
    //   - not in the set of interesting configs (=configQueue)
    currentSamples.erase(
        std::remove_if(currentSamples.begin(), currentSamples.end(),
                       [&](const auto &config) {
                         if (_invalidConfigs.count(config) > 0 or
                             std::find(configQueue.begin(), configQueue.end(), config) == configQueue.end()) {
                           return true;
                         } else {
                           return false;
                         }
                       }),
        currentSamples.end());

    // No valid configuration. This should rarely happen.
    if (currentSamples.empty()) {
      AutoPasLog(DEBUG, "Tuning could not generate a valid configuration on attempt {} of {}.", i, maxAttempts);
      if (i == maxAttempts) {
        // abort if in none of the attempts any good configuration was chosen.
        utils::ExceptionHandler::exception("BayesianSearch: Failed to sample an valid FeatureVector after {} attempts.",
                                           maxAttempts);
      }
    } else {
      // replace the config queue by what is left of the proposed configurations. Best goes to the back.
      configQueue.clear();
      std::copy(currentSamples.rbegin(), currentSamples.rend(), std::back_inserter(configQueue));
      break;
    }
  }
  return false;
}

std::vector<autopas::FeatureVector> autopas::BayesianSearch::sampleAcquisitions(size_t n,
                                                                                AcquisitionFunctionOption af) {
  // create n lhs samples
  auto currentSamples = _encoder.lhsSampleFeatures(n, _rng);

  // map container and calculate all acquisition function values
  // Use this step to also filter out any duplicate samples.
  std::map<FeatureVector, double> acquisitions;
  for (auto samplesIter = currentSamples.begin(); samplesIter != currentSamples.end();) {
    if (acquisitions.count(*samplesIter) == 0) {
      acquisitions[*samplesIter] = _gaussianProcess.calcAcquisition(af, _encoder.oneHotEncode(*samplesIter));
      ++samplesIter;
    } else {
      currentSamples.erase(samplesIter);
    }
  }

  // sort by acquisition
  std::sort(currentSamples.begin(), currentSamples.end(),
            [&acquisitions](const FeatureVector &f1, const FeatureVector &f2) {
              return acquisitions[f1] < acquisitions[f2];
            });

  return currentSamples;
}

bool autopas::BayesianSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerTraversalEstimatorOptions.empty() or
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or _dataLayoutOptions.empty() or
         _newton3Options.empty();
}

void autopas::BayesianSearch::rejectConfiguration(const autopas::Configuration &configuration, bool indefinitely) {
  if (indefinitely) {
    _invalidConfigs.insert(configuration);
  }
}

void autopas::BayesianSearch::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
  // represent a maximization problem
  _gaussianProcess.addEvidence(_encoder.oneHotEncode(configuration), -evidence.value * secondsPerMicroseconds, true);
}

bool autopas::BayesianSearch::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                    const autopas::EvidenceCollection &evidenceCollection) {
  _gaussianProcess.clear();
  return optimizeSuggestions(configQueue, evidenceCollection);
}

bool autopas::BayesianSearch::needsDomainSimilarityStatistics() const { return false; }

autopas::TuningStrategyOption autopas::BayesianSearch::getOptionType() const {
  return autopas::TuningStrategyOption::bayesianSearch;
}
