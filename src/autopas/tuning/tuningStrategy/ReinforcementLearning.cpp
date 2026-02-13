/**
 * @file ReinforcementLearning.cpp
 * @author P. Metscher
 * @date 07.05.2025
 */

#include "ReinforcementLearning.h"

#include "autopas/utils/Random.h"

autopas::ReinforcementLearning::ReinforcementLearning(const std::set<Configuration> &searchSpace,
                                                      const double learningRate, const double discountFactor,
                                                      const size_t numRandomExplorations)
    : _learningRate(learningRate), _discountFactor(discountFactor), _numRandomExplorations(numRandomExplorations) {
  if (searchSpace.size() <= _numRandomExplorations) {
    AutoPasLog(WARN, "The number of random explorations is larger than or equal to the search space size.");
  }

  if (learningRate <= 0 or learningRate > 1) {
    utils::ExceptionHandler::exception("The learning rate must be between 0 and 1 for Reinforcement Learning Tuning.");
  }

  if (discountFactor <= 0 or discountFactor > 1) {
    utils::ExceptionHandler::exception(
        "The discount factor must be between 0 and 1 for Reinforcement Learning Tuning.");
  }

  _learningRate = learningRate;
  _discountFactor = discountFactor;
}

autopas::TuningStrategyOption autopas::ReinforcementLearning::getOptionType() const {
  return TuningStrategyOption::reinforcementLearning;
}

void autopas::ReinforcementLearning::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  if (_firstTuningPhase) {
    _state.emplace(configuration, -evidence.value);
  } else {
    double oldState = _state.at(configuration);
    _state.at(configuration) = oldState + _learningRate * (-evidence.value + _discountFactor * oldState - oldState);
  }

  if (evidence.value <= _bestEvidence) {
    _bestEvidence = evidence.value;
    _bestConfiguration = configuration;
  }
}

bool autopas::ReinforcementLearning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                         const EvidenceCollection &evidenceCollection) {
  // Perform a full search over all configurations in the first tuning phase. This runs until the queue is empty.
  if (_firstTuningPhase) {
    return true;
  }

  // Perform a random search over the configurations in the exploration phase.
  if (configQueue.size() >= _numRandomExplorations and _explorationPhase) {
    // Initialize the exploration queue
    static Random randomEngine{0xF20D90B2C0FEA7AC};

    // Remove the prior best configuration from the queue
    const auto bestIt = std::ranges::find(configQueue, _bestConfiguration);
    bool hasBest = bestIt != configQueue.end();
    size_t randomSamples = _numRandomExplorations + 1;

    if (hasBest) {
      const size_t bestIdx = std::distance(configQueue.begin(), bestIt);
      std::swap(configQueue[bestIdx], configQueue[configQueue.size() - 1]);
      configQueue.resize(configQueue.size() - 1);
      randomSamples--;
    }

    // Get random elements from the queue
    std::set<Configuration> randomSearchElements =
        randomEngine.randomSubset(std::set(configQueue.begin(), configQueue.end()), randomSamples);

    configQueue.clear();
    configQueue.insert(configQueue.end(), randomSearchElements.begin(), randomSearchElements.end());

    if (hasBest) {
      configQueue.push_back(_bestConfiguration);
    }

    return false;
  }

  // Perform a greedy search over the configurations in the exploitation phase.
  if (configQueue.empty() and _explorationPhase) {
    _explorationPhase = false;

    // Find the configuration with the highest preference
    double highestEvidence = -std::numeric_limits<double>::infinity();
    Configuration preferredConfiguration;

    for (auto [config, evidence] : _state) {
      if (evidence > highestEvidence) {
        highestEvidence = evidence;
        preferredConfiguration = config;
      }
    }

    // Add the preferred configuration to the queue
    configQueue.push_back(preferredConfiguration);

    return false;
  }

  // Return true if the tuning is finished.
  if (configQueue.empty() and not _explorationPhase) {
    return true;
  }

  // Return false if the tuning is not finished, but the queue is empty.
  return false;
}

bool autopas::ReinforcementLearning::reset(size_t iteration, size_t tuningPhase,
                                           std::vector<Configuration> &configQueue,
                                           const EvidenceCollection &evidenceCollection) {
  if (configQueue.size() <= _numRandomExplorations) {
    AutoPasLog(WARN,
               "The number of random explorations is larger than or equal to the config queue size during reset.");
    return false;
  }

  _firstTuningPhase = tuningPhase == 0;
  _explorationPhase = true;
  _bestEvidence = std::numeric_limits<long>::max();

  return optimizeSuggestions(configQueue, evidenceCollection);
}
