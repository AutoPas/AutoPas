/**
 * @file ReinforcementLearning.cpp
 * @author P. Metscher
 * @date 07.05.2025
 */

#include "ReinforcementLearning.h"

#include "autopas/utils/Random.h"

autopas::ReinforcementLearning::ReinforcementLearning(const std::set<Configuration> &searchSpace,
                                                      const double learningRate, const double discountFactor,
                                                      const size_t randomExplorations)
    : _learningRate(learningRate), _discountFactor(discountFactor), _randomExplorations(randomExplorations) {
  if (searchSpace.size() <= _randomExplorations) {
    utils::ExceptionHandler::exception(
        "The search space must contain more configurations than the number of random exploration samples for "
        "Reinforcement Learning Tuning.");
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
  if (configQueue.size() >= _randomExplorations) {
    // Initialize the exploration queue
    Random randomEngine{};

    // Remove the prior best configuration from the queue
    const size_t bestIdx = std::find(configQueue.begin(), configQueue.end(), _bestConfiguration) - configQueue.begin();
    std::swap(configQueue[bestIdx], configQueue[configQueue.size() - 1]);
    configQueue.resize(configQueue.size() - 1);

    // Get random elements from the queue
    std::set<Configuration> randomSearchElements =
        randomEngine.randomSubset(std::set(configQueue.begin(), configQueue.end()), _randomExplorations);

    configQueue.clear();
    configQueue.insert(configQueue.end(), randomSearchElements.begin(), randomSearchElements.end());
    configQueue.push_back(_bestConfiguration);

    return false;
  }

  // Perform a greedy search over the configurations in the exploitation phase.
  if (configQueue.size() == 0 && _explorationPhase) {
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
  if (configQueue.size() == 0 && !_explorationPhase) {
    return true;
  }

  // Return false if the tuning is not finished, but the queue is empty.
  return false;
}

bool autopas::ReinforcementLearning::reset(size_t iteration, size_t tuningPhase,
                                           std::vector<Configuration> &configQueue,
                                           const EvidenceCollection &evidenceCollection) {
  if (configQueue.size() <= _randomExplorations) {
    utils::ExceptionHandler::exception(
        "The search space must contain more configurations than the number of random exploration samples for "
        "Reinforcement Learning Tuning.");
  }

  _firstTuningPhase = tuningPhase == 0;
  _explorationPhase = true;
  _bestEvidence = std::numeric_limits<long>::max();

  optimizeSuggestions(configQueue, evidenceCollection);

  // The first search may never contain an empty config queue.
  return false;
}
