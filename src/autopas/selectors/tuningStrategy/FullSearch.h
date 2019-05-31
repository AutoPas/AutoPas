/**
 * @file FullSearch.h
 * @author F. Gratl
 * @date 5/29/19
 */

#pragma once

#include <autopas/selectors/ContainerSelector.h>
#include <autopas/utils/ExceptionHandler.h>
#include <set>
#include "TuningStrategyInterface.h"

namespace autopas {

class FullSearch : TuningStrategyInterface {
  template <class CSelector>
  void generateSearchSpace(CSelector *containerSelector, const std::set<ContainerOption> &allowedContainerOptions,
                           std::set<TraversalOption> allowedTraversalOptions,
                           const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                           const std::set<Newton3Option> &allowedNewton3Options) {
    auto containerBefore = containerSelector->getCurrentContainer()->getContainerType();

    //@TODO needed until all containers support propper traversals
    allowedTraversalOptions.insert(TraversalOption::dummyTraversal);

    // generate all potential configs
    for (auto &containerOption : allowedContainerOptions) {
      containerSelector->selectContainer(containerOption);

      // get all traversals of the container and restrict them to the allowed ones
      auto allContainerTraversals = containerSelector->getCurrentContainer()->getAllTraversals();
      std::vector<TraversalOption> allowedAndApplicable;
      std::sort(allContainerTraversals.begin(), allContainerTraversals.end());
      std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                            allContainerTraversals.begin(), allContainerTraversals.end(),
                            std::back_inserter(allowedAndApplicable));

      for (auto &traversalOption : allowedAndApplicable) {
        for (auto &dataLayoutOption : allowedDataLayoutOptions) {
          for (auto &newton3Option : allowedNewton3Options) {
            _searchSpace.emplace(containerOption, traversalOption, dataLayoutOption, newton3Option);
          }
        }
      }
    }

    if (_searchSpace.empty()) {
      autopas::utils::ExceptionHandler::exception("FullSearch: No valid configurations could be created.");
    }

    _currentConfig = _searchSpace.begin();

    containerSelector->selectContainer(_currentConfig->_container);
  }

  Configuration getCurrentConfig() { return *_currentConfig; }

  Configuration getNextConfig() { return *(++_currentConfig); }

  std::set<Configuration> _searchSpace;
  std::set<Configuration>::iterator _currentConfig;
};
}  // namespace autopas