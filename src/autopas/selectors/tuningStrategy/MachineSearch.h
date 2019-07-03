/**
 * @file MachineSearch.h
 * @date 6/16/19
 */

#pragma once

#include <set>
#include <sstream>
#include "TuningStrategyInterface.h"
#include "autopas/AutoPas.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"
#include "fdeep/fdeep.hpp"

namespace autopas {

/**
 * Exhaustive full search of the search space by testing every applicable configuration and then selecting the optimum.
 */
template <typename Particle, typename ParticleCell>
class MachineSearch : public TuningStrategyInterface<Particle, ParticleCell> {
 public:
  /**
   * Constructor for the MachineSearch that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  MachineSearch(const std::set<ContainerOption> &allowedContainerOptions,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options, const std::string modelLink)
      : _containerOptions(allowedContainerOptions), _modelLink(modelLink) {
    // sets search space and current config
    populateSearchSpace(allowedContainerOptions, allowedTraversalOptions, allowedDataLayoutOptions,
                        allowedNewton3Options);
  }

  inline Configuration getCurrentConfiguration() override { return *_currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time) override { _traversalTimes[*_currentConfig] = time; }

  inline void reset() override {
    _traversalTimes.clear();
    _currentConfig = _searchSpace.begin();
    _configCounter = 0;
    generateMLPredictions(_modelLink);
  }

  inline bool tune() override;

  inline std::set<ContainerOption> getAllowedContainerOptions() override { return _containerOptions; };

  inline bool searchSpaceIsTrivial() override { return _searchSpace.size() == 1; }

  inline bool searchSpaceIsEmpty() override { return _searchSpace.empty(); }

  void addContainerSelector(/*const*/ ContainerSelector<Particle, ParticleCell> &containerSelector) override {
    _containerSelector = &containerSelector;
  }

 private:
  /*
   * Uses the trained ML model to choose configurations
   */
  void generateMLPredictions(std::string file);

  /*
   * Finds the next applicable ML suggestion (including current config)
   */
  void findNextSuggestion();

  /**
   * Fills the search space with the cartesian product of the given options (minus invalid combinations).
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  inline void populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                  const std::set<TraversalOption> &allowedTraversalOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options);

  inline void selectOptimalConfiguration();

  ContainerSelector<Particle, ParticleCell> *_containerSelector{nullptr};
  std::set<ContainerOption> _containerOptions;
  std::set<Configuration> _searchSpace;
  std::set<Configuration>::iterator _currentConfig;
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
  int _mlSuggestions[5];
  int _configCounter;
  double _particleCount, _boxLength, _cutoff, _verletSkin;
  bool _hasStatistics{false};
  std::string _modelLink;
};

template <typename Particle, typename ParticleCell>
void MachineSearch<Particle, ParticleCell>::findNextSuggestion() {
  while (_currentConfig != _searchSpace.end()) {
    bool isSuggestion = std::count(_mlSuggestions, _mlSuggestions + 5, _configCounter) != 0;
    if (isSuggestion) break;
    ++_currentConfig;
    ++_configCounter;
  }
}

template <typename Particle, typename ParticleCell>
void MachineSearch<Particle, ParticleCell>::generateMLPredictions(std::string file) {
  // at the start of each configuration, fill a global array with n elements for most efficient configurations
  const auto model = fdeep::load_model(file);
  double max_particle_count = 125000, max_box_length = 12, max_cutoff = 4, max_v_skin_rad = 0.3;
  _particleCount = _containerSelector->getCurrentContainer()->getNumParticles();
  _boxLength = _containerSelector->getCurrentContainer()->getBoxMax()[0] -
               _containerSelector->getCurrentContainer()->getBoxMin()[0];
  _cutoff = _containerSelector->getCurrentContainer()->getCutoff();
  _verletSkin = _containerSelector->getCurrentContainer()->getCutoff();

  const auto result = model.predict({fdeep::tensor5(
      fdeep::shape5(1, 1, 1, 1, 4),
      {static_cast<float>(_particleCount / max_particle_count), static_cast<float>(_boxLength / max_box_length),
       static_cast<float>(_cutoff / max_cutoff), static_cast<float>(_verletSkin / max_v_skin_rad)})});
  std::cout << fdeep::show_tensor5s(result) << std::endl;
  std::vector<float> probabilityVector = *result[0].as_vector();

  std::cout << "ML Suggestions are:" << std::endl;
  // find 5 largest values in propabilityVector
  for (int i = 0; i < 5; ++i) {
    // auto index = find_min(propabilityVector);
    auto index = std::max_element(probabilityVector.begin(), probabilityVector.end()) - probabilityVector.begin();
    _mlSuggestions[i] = index;
    probabilityVector[index] = 0;  // set to 0, so we don't find this again.
    std::cout << _mlSuggestions[i] << std::endl;
  }
}

template <typename Particle, typename ParticleCell>
void MachineSearch<Particle, ParticleCell>::populateSearchSpace(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options) {
  //@TODO dummyTraversal needed until all containers support propper traversals
  auto dummySet = {TraversalOption::dummyTraversal};
  std::set<TraversalOption> allowedTraversalOptionsPlusDummy;
  std::set_union(allowedTraversalOptions.begin(), allowedTraversalOptions.end(), dummySet.begin(), dummySet.end(),
                 std::inserter(allowedTraversalOptionsPlusDummy, allowedTraversalOptionsPlusDummy.begin()));

  // generate all potential configs
  for (auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    std::set<TraversalOption> allContainerTraversals = compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptionsPlusDummy.begin(), allowedTraversalOptionsPlusDummy.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (auto &traversalOption : allowedAndApplicable) {
      for (auto &dataLayoutOption : allowedDataLayoutOptions) {
        for (auto &newton3Option : allowedNewton3Options) {
          _searchSpace.emplace(containerOption, traversalOption, dataLayoutOption, newton3Option);
        }
      }
    }
  }

  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("MachineSearch: No valid configurations could be created.");
  }

  _currentConfig = _searchSpace.begin();
  _configCounter = 0;
  // generateMLPredictions(_modelLink);
  // findNextSuggestion();
}

template <typename Particle, typename ParticleCell>
bool MachineSearch<Particle, ParticleCell>::tune() {
  AutoPasLog(debug, "You are in MachineSearch::tune()");
  // repeat as long as traversals are not applicable or we run out of configs
  ++_currentConfig;
  ++_configCounter;
  findNextSuggestion();
  if (_currentConfig == _searchSpace.end()) {
    selectOptimalConfiguration();
    return false;
  }
  return true;
}

template <typename Particle, typename ParticleCell>
void MachineSearch<Particle, ParticleCell>::selectOptimalConfiguration() {
  if (_searchSpace.size() == 1) {
    _currentConfig = _searchSpace.begin();
    return;
  }

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "MachineSearch: Trying to determine fastest configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  auto optimum = std::min_element(_traversalTimes.begin(), _traversalTimes.end(),
                                  [](std::pair<Configuration, size_t> a, std::pair<Configuration, size_t> b) -> bool {
                                    return a.second < b.second;
                                  });

  _currentConfig = _searchSpace.find(optimum->first);
  // sanity check
  if (_currentConfig == _searchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "MachineSearch: Optimal configuration not found in list of configurations!");
  }

  // measurements are not needed anymore
  _traversalTimes.clear();

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}

template <typename Particle, typename ParticleCell>
void MachineSearch<Particle, ParticleCell>::removeN3Option(Newton3Option badNewton3Option) {
  for (auto ssIter = _searchSpace.begin(); ssIter != _searchSpace.end();) {
    if (ssIter->_newton3 == badNewton3Option) {
      // change current config to the next non-deleted
      if (ssIter == _currentConfig) {
        ssIter = _searchSpace.erase(ssIter);
        _currentConfig = ssIter;
      } else {
        ssIter = _searchSpace.erase(ssIter);
      }
    } else {
      ++ssIter;
    }
  }

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

}  // namespace autopas
