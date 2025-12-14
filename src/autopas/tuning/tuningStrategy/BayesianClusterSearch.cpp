/**
 * @file BayesianClusterSearch.cpp
 * @author Fabio Gratl
 * @date 17.11.22
 */

#include "BayesianClusterSearch.h"

#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"

autopas::BayesianClusterSearch::BayesianClusterSearch(
    const InteractionTypeOption &interactionType, const std::set<ContainerOption> &allowedContainerOptions,
    const NumberSet<double> &allowedCellSizeFactors, const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options,
    const std::set<VectorizationPatternOption> &allowedVecPatternOptions, size_t maxEvidence,
    AcquisitionFunctionOption predAcqFunction, const std::string &outputSuffix, size_t predNumLHSamples,
    unsigned long seed)
    : _interactionType(interactionType),
      _containerOptionsSet(allowedContainerOptions),
      _vecPatternOptions(allowedVecPatternOptions.begin(), allowedVecPatternOptions.end()),
      _dataLayoutOptions(allowedDataLayoutOptions.begin(), allowedDataLayoutOptions.end()),
      _newton3Options(allowedNewton3Options.begin(), allowedNewton3Options.end()),
      _cellSizeFactors(allowedCellSizeFactors.clone()),
      _encoder(),
      _invalidConfigs(),
      _rng(seed),
      _gaussianCluster({}, continuousDims, GaussianCluster::WeightFunction::evidenceMatchingScaledProbabilityGM, sigma,
                       _rng, GaussianCluster::defaultVecToString, outputSuffix),
      _neighbourFun([this](const Eigen::VectorXi &target) -> std::vector<std::pair<Eigen::VectorXi, double>> {
        return _encoder.clusterNeighboursManhattan1Container(target);
      }),
      _maxEvidence(maxEvidence),
      _predAcqFunction(predAcqFunction),
      _predNumLHSamples(predNumLHSamples),
      _firstTuningPhase(true),
      _currentIteration(0),
      _iterationScale(1. / (maxEvidence * iterationScalePerMaxEvidence)),
      _currentNumEvidence(0),
      _currentOptimalTime(std::numeric_limits<long>::max()) {
  if (predNumLHSamples <= 0) {
    utils::ExceptionHandler::exception(
        "BayesianClusterSearch: Number of samples used for predictions must be greater than 0!");
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
    autopas::utils::ExceptionHandler::exception("BayesianClusterSearch: No valid configurations could be created.");
  }

  updateOptions();
  _gaussianCluster.setVectorToStringFun(
      [this](const GaussianModelTypes::VectorPairDiscreteContinuous &vec) -> std::string {
        return _encoder.convertFromCluster(vec).toString();
      });
}

autopas::BayesianClusterSearch::~BayesianClusterSearch() = default;

void autopas::BayesianClusterSearch::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  // store optimal evidence
  if (evidence.value < _currentOptimalTime) {
    _currentOptimalTime = evidence.value;
    _currentOptimalConfig = configuration;
  }

  // encoded vector
  const auto vec = _encoder.convertToCluster(configuration, evidence.iteration * _iterationScale);
  // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
  // represent a maximization problem
  _gaussianCluster.addEvidence(vec, -evidence.value * secondsPerMicroseconds);

  _currentIteration = evidence.iteration;
  ++_currentNumEvidence;
  _currentAcquisitions.clear();
}

bool autopas::BayesianClusterSearch::reset(size_t iteration, size_t tuningPhase,
                                           std::vector<Configuration> &configQueue,
                                           const autopas::EvidenceCollection &evidenceCollection) {
  const auto iterationSinceLastEvidence = iteration - _currentIteration;
  if (static_cast<double>(iterationSinceLastEvidence) * _iterationScale > suggestedMaxDistance) {
    AutoPasLog(WARN,
               "BayesianClusterSearch: Time since the last evidence may be too long ({} > {}). "
               "You should either decrease the number of iterations between tuning phases "
               "or gather more evidence (currently: {}).",
               iterationSinceLastEvidence, suggestedMaxDistance / _iterationScale, _maxEvidence);
  }

  _currentIteration = iteration;
  _currentNumEvidence = 0;
  _currentOptimalTime = std::numeric_limits<long>::max();
  _currentAcquisitions.clear();

  _firstTuningPhase = tuningPhase == 0 ? true : false;

  if (not _firstTuningPhase) {
    return optimizeSuggestions(configQueue, evidenceCollection);
  }

  // BayesianClusterSearch does no intentional config wipes to stop the tuning phase
  return false;
}

bool autopas::BayesianClusterSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerTraversalEstimatorOptions.empty() or
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or _dataLayoutOptions.empty() or
         _newton3Options.empty();
}

void autopas::BayesianClusterSearch::updateOptions() {
  _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options,
                             *_cellSizeFactors, _vecPatternOptions);

  auto newRestrictions = _encoder.getDiscreteRestrictions();
  _gaussianCluster.setDimensions(std::vector<int>(newRestrictions.begin(), newRestrictions.end()));
}

bool autopas::BayesianClusterSearch::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                         const EvidenceCollection &evidenceCollection) {
  // In the first tuning phase do nothing since we first need some data.
  if (_firstTuningPhase) {
    return false;
  }

  // no more tunings steps
  if (_currentNumEvidence >= _maxEvidence) {
    // select best config of current tuning phase
    configQueue.clear();
    return true;
  }

  // try to sample a valid vector which is expected to yield a good acquisition
  for (size_t i = 0; i < maxAttempts; ++i) {
    auto currentAcquisitions = sampleAcquisitions(_predNumLHSamples, _predAcqFunction);

    // Filter out all acquisitions, which do map to a configuration which is either:
    //   - invalid
    //   - not in the set of interesting configs (=configQueue)
    currentAcquisitions.erase(
        std::remove_if(currentAcquisitions.begin(), currentAcquisitions.end(),
                       [&](const auto &acquisition) {
                         const auto config = _encoder.convertFromCluster(acquisition);
                         if (_invalidConfigs.count(config) > 0 or
                             std::find(configQueue.begin(), configQueue.end(), config) == configQueue.end()) {
                           return true;
                         } else {
                           return false;
                         }
                       }),
        currentAcquisitions.end());

    // No valid configuration. This should rarely happen.
    if (currentAcquisitions.empty()) {
      AutoPasLog(DEBUG, "Tuning could not generate a valid configuration on attempt {} of {}.", i, maxAttempts);
      if (i == maxAttempts) {
        utils::ExceptionHandler::exception(
            "BayesianClusterSearch: Failed to sample an valid FeatureVector after {} attempts.", maxAttempts);
      }
    } else {
      // replace the config queue by what is left of the proposed configurations.
      configQueue.clear();
      std::transform(currentAcquisitions.rbegin(), currentAcquisitions.rend(), std::back_inserter(configQueue),
                     [&](const auto &acquisition) { return _encoder.convertFromCluster(acquisition); });
      break;
    }
  }
  return false;
}

std::vector<autopas::GaussianModelTypes::VectorPairDiscreteContinuous>
autopas::BayesianClusterSearch::sampleAcquisitions(size_t n, AcquisitionFunctionOption af) {
  // create n lhs samples
  auto continuousSamples = _encoder.lhsSampleFeatureCluster(n, _rng, _currentIteration * _iterationScale);

  // calculate all acquisitions
  const auto currentAcquisitions = _gaussianCluster.sampleOrderedByAcquisition(af, _neighbourFun, continuousSamples);
  return currentAcquisitions;
}

void autopas::BayesianClusterSearch::rejectConfiguration(const autopas::Configuration &configuration,
                                                         bool indefinitely) {
  _invalidConfigs.insert(configuration);
}
autopas::TuningStrategyOption autopas::BayesianClusterSearch::getOptionType() const {
  return TuningStrategyOption::bayesianClusterSearch;
}
