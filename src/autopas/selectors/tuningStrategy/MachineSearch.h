/**
 * @file MachineSearch.h
 * @author candas
 * @date 6/16/19
 */

#pragma once
#include <fdeep/fdeep.hpp>
#include <set>
#include <sstream>
#include "TuningStrategyInterface.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"
#include <fstream>
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * Exhaustive full search of the search space by testing every applicable configuration and then selecting the optimum.
 */
template <class Particle, class ParticleCell>
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
                const std::set<double> &allowedCellSizeFactors,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options, const std::string& modelLink)
      : _containerOptions(allowedContainerOptions){
    // sets search space and current config
    populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                        allowedDataLayoutOptions, allowedNewton3Options);

    // Here we read our configurations for our model
    std::ifstream confFile;
    confFile.open(modelLink);
    if (!confFile) {
        std::cout << "Unable to open file " << modelLink << std::endl;
        exit(1);
    }
    std::string line;
    int lineNum = 0;
    while (std::getline(confFile, line)) {
      if (line.empty() || line.at(0) == '#') { continue; }
      lineNum++;
      if (lineNum == 1) {_mlmodel = std::make_unique<fdeep::model>(fdeep::load_model(path));}
      else if (lineNum == 2) {_inputConfig = std::stoi(line);}
      else if (lineNum == 3 && line.substr(0, 5) == "Norm:") {
        _normalize = true;
        auto inputs = autopas::utils::StringUtils::tokenize(line.substr(5, line.size() - 5), " ,;|/");
        switch(_inputConfig) {
          case 2 : _picture = std::stoi(inputs[0]);
                   break;
          case 3 :
          case 1 : _particleCount = std::stod(inputs[0]);
                   _boxLength = std::stod(inputs[1]);
                   _cutoff = std::stod(inputs[2]);
                   _verletSkin = std::stod(inputs[3]);
                   if (_inputConfig == 3) { _picture = std::stoi(inputs[4]); }
        }
      }
      else {
        auto outputs = autopas::utils::StringUtils::tokenize(line, ",;|/"); // don't tokenize with whitespace
        ContainerOption container = *autopas::utils::StringUtils::parseContainerOptions(outputs[1]).begin();
        TraversalOption traversal = *autopas::utils::StringUtils::parseTraversalOptions(outputs[3]).begin();
        DataLayoutOption dataLayout = *autopas::utils::StringUtils::parseDataLayout(outputs[5]).begin();
        Newton3Option newton3 = *autopas::utils::StringUtils::parseNewton3Options(outputs[7]).begin();
        _outputConfig.push_back(std::make_tuple(container, 1., traversal, dataLayout, newton3));
      }
    }

    // TODO: Maybe also check if all option sets really have only one element
    /* Iterating through the configurations would be nice to check if everything went smoothly
    for (auto ocIter = _outputConfig.begin(); ocIter != _outputConfig.end();) {
      std::cout << ocIter->toString() << std::endl;
      ++ocIter;
    } */

    // Error checking
    if (1 > _inputConfig || _inputConfig > 3) {
        std::cout << "Input configuration must be between 1 and 3" << std::endl;
        exit(1);
    }

    if (_normalize && (_inputConfig == 1 || _inputConfig == 3) &&
    (_particleCount <= 0 || _boxLength <= 0 || _cutoff <= 0 || _verletSkin <= 0)) {
        std::cout << "Normalization was requested with faulty values" << std::endl;
        exit(1);
    }

    if (_normalize && (_inputConfig == 2 || _inputConfig == 3) && (_picture <= 0)) {
        std::cout << "Normalization was requested with faulty values" << std::endl;
        exit(1);
    }

    //_mlmodel = fdeep::load_model("fake");

    // Generate an output vector to make sure configuration makes sense
    const auto result = _mlmodel->predict({fdeep::tensor5(fdeep::shape5(1, 1, 1, 1, 4), {0, 0, 0, 0})});
    std::vector<float> probabilityVector = *result[0].as_vector();

    if (_outputConfig.size() != probabilityVector.size()) {
        std::cout << "Output configuration must have the same size as the model's output" << std::endl;
        exit(1);
    }

  }

  inline const Configuration &getCurrentConfiguration() override { return *_currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time) override { _traversalTimes[*_currentConfig] = time; }

  inline void reset() override {
    _traversalTimes.clear();
    //_currentConfig = _searchSpace.end();
    //_configCounter = 0;
    generateMLPredictions();
    _currentConfig = _searchSpace.begin();
    for (auto ssIter = _searchSpace.begin(); ssIter != _searchSpace.end();) {
        std::cout << ssIter->toString() << std::endl;
        ++ssIter;
    }

  }

  inline bool tune(bool = false) override;

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
  void generateMLPredictions();

  /*
   * Finds the next applicable ML suggestion (including current config)
   */
  //void findNextSuggestion();

  /**
   * Fills the search space with the cartesian product of the given options (minus invalid combinations).
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  inline void populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                  const std::set<double> &allowedCellSizeFactors,
                                  const std::set<TraversalOption> &allowedTraversalOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options);

  inline void selectOptimalConfiguration();

  std::set<ContainerOption> _containerOptions;
  std::set<Configuration> _searchSpace;
  std::set<Configuration>::iterator _currentConfig;
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
  ContainerSelector<Particle, ParticleCell> *_containerSelector;

  //int _configCounter;
  std::array<int, 5> _mlSuggestions;
  // the following doubles and _picture are maximum values to normalize with
  double _particleCount, _boxLength, _cutoff, _verletSkin; // assume that they are initialized with 0.0
  int _picture = 0;
  std::unique_ptr<fdeep::model> _mlmodel{nullptr};
  int _inputConfig = 0;
  std::vector<std::tuple<ContainerOption, double, TraversalOption, DataLayoutOption, Newton3Option>> _outputConfig;
  bool _normalize = false;
};

/*
template <class Particle, class ParticleCell>
void MachineSearch<Particle, ParticleCell>::findNextSuggestion() {
  if (_currentConfig == _searchSpace.end()) {
    _currentConfig = _searchSpace.begin();
  } else {
    ++_currentConfig;
    ++_configCounter;
  }
  while (_currentConfig != _searchSpace.end()) {
    bool isSuggestion = std::count(_mlSuggestions, _mlSuggestions + 5, _configCounter) != 0;
    if (isSuggestion) break;
    ++_currentConfig;
    ++_configCounter;
  }
}*/

template <class Particle, class ParticleCell>
void MachineSearch<Particle, ParticleCell>::generateMLPredictions() {
  // at the start of each configuration, fill a global array with n elements for most efficient configurations
  double max_particle_count = 125000, max_box_length = 12, max_cutoff = 4, max_v_skin_rad = 0.3;
  auto container = _containerSelector->getCurrentContainer();
  _particleCount = container->getNumParticles();
  _boxLength = container->getBoxMax()[0] - container->getBoxMin()[0];
  _cutoff = container->getCutoff();
  _verletSkin = container->getSkin();

  float in_par = (_particleCount > max_particle_count) ? 1 : _particleCount / max_particle_count;
  float in_box = (_boxLength > max_box_length) ? 1 : _boxLength / max_box_length;
  float in_cut = (_cutoff > max_cutoff) ? 1 : _cutoff / max_cutoff;
  float in_ver = (_verletSkin > max_v_skin_rad) ? 1 : _verletSkin / max_v_skin_rad;

  const auto result =
      _mlmodel->predict({fdeep::tensor5(fdeep::shape5(1, 1, 1, 1, 4), {in_par, in_box, in_cut, in_ver})});
  std::cout << fdeep::show_tensor5s(result) << std::endl;
  std::vector<float> probabilityVector = *result[0].as_vector();

  std::cout << "ML Suggestions are:" << std::endl;
  // find 5 largest values in propabilityVector
  for (int i = 0; i < 5; ++i) {
    // auto index = find_min(propabilityVector);
    auto index = std::max_element(probabilityVector.begin(), probabilityVector.end()) - probabilityVector.begin();
    _mlSuggestions[i] = index;
    probabilityVector[index] = 0;  // set to 0, so we don't find this again.
    std::cout << _mlSuggestions[i] << ", ";
  }
  std::cout << std::endl;
  _searchSpace.clear();  // empty the set

  for (int i = 0; i < 5; ++i) {
    auto tuple = _outputConfig[_mlSuggestions[i]];
    _searchSpace.emplace(std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple), std::get<3>(tuple), std::get<4>(tuple));
    // std::apply(_searchSpace.emplace, tuple); // is much more readable, but error
  }
}

template <class Particle, class ParticleCell>
void MachineSearch<Particle, ParticleCell>::populateSearchSpace(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
    const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options) {
    // generate all potential configs
    for (auto &containerOption : allowedContainerOptions) {
        // get all traversals of the container and restrict them to the allowed ones
        const std::set<TraversalOption> &allContainerTraversals =
                compatibleTraversals::allCompatibleTraversals(containerOption);
        std::set<TraversalOption> allowedAndApplicable;
        std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                              allContainerTraversals.begin(), allContainerTraversals.end(),
                              std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

        for (auto &cellSizeFactor : allowedCellSizeFactors)
            for (auto &traversalOption : allowedAndApplicable) {
                for (auto &dataLayoutOption : allowedDataLayoutOptions) {
                    for (auto &newton3Option : allowedNewton3Options) {
                        _searchSpace.emplace(containerOption, cellSizeFactor, traversalOption, dataLayoutOption, newton3Option);
                    }
                }
            }
    }

    AutoPasLog(debug, "Points in search space: {}", _searchSpace.size());

  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("MachineSearch: No valid configurations could be created.");
  }

  _currentConfig = _searchSpace.begin();  // this has to be valid, to add particles to the container
  //_configCounter = 0;
  // generateMLPredictions(_modelLink);
  // findNextSuggestion();
}

template <class Particle, class ParticleCell>
bool MachineSearch<Particle, ParticleCell>::tune(bool) {
  AutoPasLog(debug, "You are in MachineSearch::tune()");
  // repeat as long as traversals are not applicable or we run out of configs
  //findNextSuggestion();
  ++_currentConfig;
  if (_currentConfig == _searchSpace.end()) {
    selectOptimalConfiguration();
    return false;
  }
  return true;
}

template <class Particle, class ParticleCell>
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

template <class Particle, class ParticleCell>
void MachineSearch<Particle, ParticleCell>::removeN3Option(Newton3Option badNewton3Option) {
  for (auto ssIter = _searchSpace.begin(); ssIter != _searchSpace.end();) {
    if (ssIter->newton3 == badNewton3Option) {
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
