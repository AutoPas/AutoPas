/**
 * @file SimulatedAnnealing.cpp
 * @author Niklas Ladurner
 * @date 08.11.2026
 */

#include "SimulatedAnnealing.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <optional>
#include <vector>

#include "autopas/tuning/Configuration.h"

autopas::SimulatedAnnealing::SimulatedAnnealing(size_t maxEvidence, unsigned long seed, double initialTemperature,
                                                double endTemperature, double exclusionFactor,
                                                double levelChangeVariance)
    : _maxEvidence(maxEvidence),
      _numEvidenceCollected(0),
      _rng(seed),
      _initialTemperature(initialTemperature),
      _currentTemperature(0.0),
      _endTemperature(endTemperature),
      _exclusionFactor(exclusionFactor),
      _levelChangeVariance(levelChangeVariance) {
  if (endTemperature >= initialTemperature) {
    utils::ExceptionHandler::exception(
        "SimulatedAnnealing::endTemperature has to be smaller than SimulatedAnnealing::initialTemperature ({}) but is "
        "{}.",
        initialTemperature, endTemperature);
  } else if (endTemperature <= 0) {
    utils::ExceptionHandler::exception("SimulatedAnnealing::endTemperature has to be bigger than 0.0 but is {}.",
                                       endTemperature);
  }

  if (levelChangeVariance <= 0) {
    utils::ExceptionHandler::exception("SimulatedAnnealing::levelChangeVariance has to be bigger than 0.0 but is {}.",
                                       levelChangeVariance);
  }

  if (exclusionFactor <= 0) {
    utils::ExceptionHandler::exception("SimulatedAnnealing::exclusionFactor has to be bigger than 0.0 but is {}.",
                                       exclusionFactor);
  }
}

namespace {
/// @brief Computes a "level" difference value for two configurations. Parameters that indicate larger "difference"
/// between two configurations receive a larger value. E.g. if two configurations differ in the container they use, they
/// receive the highest difference level of 3, as switching from one configuration to the other would be costly
/// (rebuilding of container). The levels range from 0 (smallest) to 3 (highest) difference.
/// @param a Configuration to compare.
/// @param b Configuration to compare.
/// @return The highest difference level of the two configurations: 3 (container), 2 (traversal), 1 (data
/// layout, newton3) or 0 (cellSizeFactor, loadEstimator, vecPattern)
int levelDiff(const autopas::Configuration &a, const autopas::Configuration &b) {
  if (a.container != b.container) return 3;
  if (a.traversal != b.traversal) return 2;
  if (a.dataLayout != b.dataLayout || a.newton3 != b.newton3) return 1;
  return 0;
}

/// @brief Compares two configurations and returns whether both configurations have the same value(s) at the given
/// difference level. Level 0 always returns true.
/// @param a Configuration to compare.
/// @param b Configuration to compare.
/// @param level Difference level to compare at.
/// @return Returns whether the configurations have the same value(s) or not.
bool sameValueAtLevel(const autopas::Configuration &a, const autopas::Configuration &b, int level) {
  if (level == 3) return a.container == b.container;
  if (level == 2) return a.traversal == b.traversal;
  if (level == 1) return (a.dataLayout == b.dataLayout) && (a.newton3 == b.newton3);
  return true;
}

/// @brief Computes a simple heuristic of the potential cost of switching from one configuration to the other, up to a
/// given difference level.
/// @param a Configuration to compare.
/// @param b Configuration to compare.
/// @param level Difference level to compare up to.
/// @return The difference heuristic between the two configurations.
int parameterDiff(const autopas::Configuration &a, const autopas::Configuration &b, int level) {
  // @todo Possibly penalize differently depending on the actual parameters.
  int differenceCounter = 0;

  differenceCounter += a.cellSizeFactor != b.cellSizeFactor ? 1 : 0;
  differenceCounter += a.loadEstimator != b.loadEstimator ? 1 : 0;
  differenceCounter += a.vecPattern != b.vecPattern ? 1 : 0;

  if (level > 1) {
    differenceCounter += a.dataLayout != b.dataLayout ? 1 : 0;
    differenceCounter += a.newton3 != b.newton3 ? 1 : 0;
  }

  if (level > 2) {
    differenceCounter += a.traversal != b.traversal ? 1 : 0;
  }

  return differenceCounter;
}

/// @brief Computes a prediction for Δ = f(proposed)/f(current) - 1 based on measurements of previous tuning phases.
/// This allows to sort out candidates that are expected to perform very poorly, so that no tuning iterations are wasted
/// on sampling candidates that won't be selected.
/// @param evidenceCollection Evidencecollection passed by the AutoTuner.
/// @param current Current configuration.
/// @param candidate Candidate configuration.
/// @param tuningPhase Current tuning phase.
/// @return The predicted value for Δ or std::nullopt if there has been no tuning phase in the lookback window in which
/// both configurations were sampled.
std::optional<double> predictDelta(const autopas::EvidenceCollection &evidenceCollection,
                                   const autopas::Configuration &current, const autopas::Configuration &candidate,
                                   const size_t tuningPhase) {
  // Maximum number of tuning phases to look back for old measurements.
  // @todo Possibly make this configurable?
  const size_t maxMeasurementAge = 2;
  const auto evidenceOfCurrentPtr = evidenceCollection.getEvidence(current);
  const auto evidenceOfCandidatePtr = evidenceCollection.getEvidence(candidate);

  if (evidenceOfCurrentPtr == nullptr || evidenceOfCandidatePtr == nullptr) return std::nullopt;

  const auto &evidenceOfCurrent = *evidenceOfCurrentPtr;
  const auto &evidenceOfCandidate = *evidenceOfCandidatePtr;

  for (size_t revPhaseIndex = 1; revPhaseIndex <= maxMeasurementAge; revPhaseIndex++) {
    if (tuningPhase < revPhaseIndex) continue;
    const size_t phaseToLookAt = tuningPhase - revPhaseIndex;
    const auto measurementOfCurrentIt =
        std::find_if(evidenceOfCurrent.begin(), evidenceOfCurrent.end(),
                     [phaseToLookAt](const autopas::Evidence &e) { return e.tuningPhase == phaseToLookAt; });
    const auto measurementOfCandidateIt =
        std::find_if(evidenceOfCandidate.begin(), evidenceOfCandidate.end(),
                     [phaseToLookAt](const autopas::Evidence &e) { return e.tuningPhase == phaseToLookAt; });

    if (measurementOfCurrentIt == evidenceOfCurrent.end() || measurementOfCandidateIt == evidenceOfCandidate.end())
      continue;

    return std::make_optional(static_cast<double>(measurementOfCandidateIt->effectiveValue) /
                                  static_cast<double>(measurementOfCurrentIt->effectiveValue) -
                              1);
  }
  return std::nullopt;
}
}  // namespace

void autopas::SimulatedAnnealing::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  ++_numEvidenceCollected;

  if (!_currentState.has_value()) {
    _currentState = std::make_optional(AlgorithmState{configuration, evidence.effectiveValue});
  } else if (_pendingProposal.has_value() && _pendingProposal.value() == configuration) {
    // Only compute delta if the measured configuration was actually our proposal.
    const double delta =
        static_cast<double>(evidence.effectiveValue) / static_cast<double>(_currentState->effectiveValue) - 1;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    const double rand = distribution(_rng);

    // Check the metropolis criterion
    const bool accept = (delta <= 0) || rand < std::exp(-delta / _currentTemperature);
    if (accept) _currentState = std::make_optional(AlgorithmState{configuration, evidence.effectiveValue});

    _currentTemperature *= _coolingFactor;
  }
}

bool autopas::SimulatedAnnealing::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                      const EvidenceCollection &evidenceCollection) {
  if (_numEvidenceCollected == _maxEvidence) {
    configQueue.clear();
    return true;
  }

  if (configQueue.empty() || !_currentState.has_value()) {
    return false;
  }

  // We select a new proposal as follows: We categorize the configuration parameter space into 4 distinct levels:
  // 3 (container), 2 (traversal), 1 (data layout, newton3) or 0 (cellSizeFactor, loadEstimator, vecPattern and other
  // fine-tuning parameters). Over the course of the tuning phase, as the temperature cools down, we use shifted
  // normal distributions to decide at which level the next proposed config should differ from the current one. With
  // that, large changes (such as changing the container) are more likely at higher temperature (i.e. in the beginning
  // of the tuning pahse) and small fine-tuning changes are more likely later in the tuning phase.

  // For that, we first sort the candidate configurations by the highest level they differ from the current
  // configuration.

  std::vector<std::vector<Configuration>> configsOrderedByLevel(4);

  for (const auto &candidate : configQueue) {
    if (candidate == _currentState->config) continue;

    // Compute a prediction of Δ on old measurements to eliminate candidates that may take up valuable evaluations.
    const auto prediction = predictDelta(evidenceCollection, _currentState->config, candidate, _currentTuningPhase);
    if (prediction.has_value() && prediction.value() > _exclusionFactor * _currentTemperature) continue;

    configsOrderedByLevel.at(levelDiff(_currentState->config, candidate)).push_back(candidate);
  }

  // Normalize the current temperature to [0, 3]. Note that the temperature falls geometrically, but we want a linear
  // progress value.
  const double normalizedProgress = std::clamp(
      3.0 * std::log(_currentTemperature / _endTemperature) / std::log(_initialTemperature / _endTemperature), 0.0,
      3.0);

  // Calculate weights (≈ probabilities) of change happening on a given level (unnormalized).
  std::vector<double> levelWeights(4);
  bool candidatesExist = false;
  for (size_t levelIndex = 0; levelIndex < levelWeights.size(); levelIndex++) {
    // If there exists no possible candidate at level i, set the probability to 0.
    // If one exists, compute the probability as a normal distribution centered at the levelIndex.
    if (configsOrderedByLevel.at(levelIndex).empty()) {
      levelWeights.at(levelIndex) = 0.0;
    } else {
      levelWeights.at(levelIndex) = std::exp(
          -std::pow((normalizedProgress - static_cast<double>(levelIndex)), 2.0) / (2.0 * _levelChangeVariance));
      candidatesExist = true;
    }
  }

  if (!candidatesExist) {
    // No promising candidates left to propose, clear queue.
    configQueue.clear();
    return true;
  }

  // Select level at which to select next proposal.
  std::discrete_distribution<int> levelDistribution(levelWeights.begin(), levelWeights.end());
  const int selectedLevelIndex = levelDistribution(_rng);

  // Select differing parameter value evenly among available values. For that, we select some configurations as
  // representatives for their value of the given parameter, as to not skew the random pick by e.g. containers which
  // have many traversal options.
  std::vector<Configuration> representatives;
  const auto &members = configsOrderedByLevel.at(selectedLevelIndex);
  std::copy_if(members.begin(), members.end(), std::back_inserter(representatives), [&](const auto &c) {
    return std::none_of(representatives.begin(), representatives.end(),
                        [&](const auto &m) { return sameValueAtLevel(c, m, selectedLevelIndex); });
  });

  // Now select uniformly among the representatives about which value should be chosen.
  std::uniform_int_distribution<std::mt19937::result_type> representativeDistribution(0, representatives.size() - 1);
  const Configuration selectedRepresentative = representatives.at(representativeDistribution(_rng));

  // Retain only candidates with the selected parameter value that also require a minimal change from the current
  // configuration in the other parameters.
  std::vector<Configuration> remainingCandidates;
  int minimalDifference = std::numeric_limits<int>::max();
  for (const auto &c : members) {
    if (sameValueAtLevel(c, selectedRepresentative, selectedLevelIndex)) {
      const int candidateDifference = parameterDiff(c, _currentState->config, selectedLevelIndex);
      if (candidateDifference < minimalDifference) {
        minimalDifference = candidateDifference;
        remainingCandidates.clear();
      }

      if (candidateDifference == minimalDifference) remainingCandidates.push_back(c);
    }
  }

  // Lastly, select uniformly among the remaining candidates
  std::uniform_int_distribution<std::mt19937::result_type> remainingDistribution(0, remainingCandidates.size() - 1);
  const auto &proposal = remainingCandidates.at(remainingDistribution(_rng));
  std::swap(configQueue.back(), *std::find(configQueue.begin(), configQueue.end(), proposal));
  _pendingProposal = std::make_optional(configQueue.back());

  return false;
}

bool autopas::SimulatedAnnealing::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                        const autopas::EvidenceCollection &evidenceCollection) {
  _currentTemperature = _initialTemperature;
  _numEvidenceCollected = 0;
  _currentState.reset();
  _pendingProposal.reset();

  _currentTuningPhase = tuningPhase;

  const size_t numEvidenceInPhase = std::min(_maxEvidence, configQueue.size());
  _coolingFactor = numEvidenceInPhase > 2 ? std::pow(_endTemperature / _initialTemperature,
                                                     1.0 / static_cast<double>(numEvidenceInPhase - 2))
                                          : 1.0;

  if (!evidenceCollection.empty()) {
    // If we already collected evidence, move the latest optimal configuration to the back
    const auto it = std::find(configQueue.begin(), configQueue.end(),
                              std::get<0>(evidenceCollection.getLatestOptimalConfiguration()));

    if (it == configQueue.end()) return false;
    std::swap(configQueue.back(), *it);
  }

  return false;
}

autopas::TuningStrategyOption autopas::SimulatedAnnealing::getOptionType() const {
  return TuningStrategyOption::simulatedAnnealing;
}
