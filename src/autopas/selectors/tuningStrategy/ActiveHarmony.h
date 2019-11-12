/**
 * @file ActiveHarmony.h
 * @author Jakob Englhauser
 * @date 21.10.19
 */

#pragma once

#include <string>
#include <autopas/selectors/Configuration.h>
#include <autopas/containers/CompatibleTraversals.h>
#include "TuningStrategyInterface.h"
#include "include/hclient.h"

namespace autopas {

class ActiveHarmony : public TuningStrategyInterface {

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  ActiveHarmony(const std::set<ContainerOption> &allowedContainerOptions,
                const NumberInterval<double> &allowedCellSizeFactors,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options)
          : _allowedContainerOptions(),
            _allowedCellSizeFactors(allowedCellSizeFactors.clone()),
            _allowedTraversalOptions(allowedTraversalOptions),
            _allowedDataLayoutOptions(allowedDataLayoutOptions),
            _allowedNewton3Options(allowedNewton3Options),
            _currentConfig() {

    // reduce traversal and container option to possible combinations
    for (auto &traversalOption : _allowedTraversalOptions) {
      auto container = *compatibleTraversals::allCompatibleContainers(traversalOption).begin();
      if (allowedContainerOptions.find(container) == allowedContainerOptions.end()) {
        _allowedTraversalOptions.erase(traversalOption);
      } else {
        _allowedContainerOptions.emplace(container);
      }
    }
    if (getenv("HARMONY_HOME") == nullptr) {
      putenv(HARMONY_HOME);
    }

    reset();
  }

  ~ActiveHarmony() override {
    if (htask != nullptr) {
      ah_leave(htask);
      ah_kill(htask);
    }
    if (hdesc != nullptr) {
      ah_close(hdesc);
      ah_free(hdesc);
    }
  }

  inline void addEvidence(long time) override;

  inline bool tune(bool currentInvalid) override;

  inline void removeN3Option(Newton3Option option) override;

  inline bool searchSpaceIsTrivial() const override;

  inline bool searchSpaceIsEmpty() const override;

  inline const Configuration &getCurrentConfiguration() const override;

  inline void reset() override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override;

 private:
  hdesc_t *hdesc;
  htask_t *htask;

  std::set<ContainerOption> _allowedContainerOptions;
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors; // Maybe use min, max, stepsize
  std::set<TraversalOption> _allowedTraversalOptions;
  std::set<DataLayoutOption> _allowedDataLayoutOptions;
  std::set<Newton3Option> _allowedNewton3Options;

  Configuration _currentConfig;

  static constexpr int cellSizeSamples = 100;
};

void ActiveHarmony::addEvidence(long time) {
  auto perf = (double) time;
  if (ah_report(htask, &perf) != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::addEvidence: Error reporting performance to server");
  }
  AutoPasLog(debug, "ActiveHarmony::addEvidence: Reported time {} for configuration {}.", perf,
             _currentConfig.toString());
}

bool ActiveHarmony::tune(bool currentInvalid) {
  if (currentInvalid) {
    auto perf = std::numeric_limits<double>::max();
    ah_report(htask, &perf);
  }
  if (ah_fetch(htask) < 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::tune: Error fetching values from server");
  }
  std::string traversalOptionStr = ah_get_enum(htask, "traversalOption");
  auto traversalOption = TraversalOption::parseOptionExact(traversalOptionStr);
  std::string dataLayoutOption = ah_get_enum(htask, "dataLayoutOption");
  std::string newton3Option = ah_get_enum(htask, "newton3Option");
  double cellSizeFactor = ah_get_real(htask, "cellSizeFactor");
  _currentConfig = Configuration(*compatibleTraversals::allCompatibleContainers(traversalOption).begin(),
                                 cellSizeFactor, traversalOption, DataLayoutOption::parseOptionExact(dataLayoutOption),
                                 Newton3Option::parseOptionExact(newton3Option));
  if (_allowedNewton3Options.find(_currentConfig.newton3) == _allowedNewton3Options.end()) {
    return tune(true);
  }
  AutoPasLog(debug, "ActiveHarmony::tune: Trying configuration {}.", _currentConfig.toString());
  auto converged = ah_converged(htask);
  if (converged) {
    AutoPasLog(debug, "ActiveHarmony::tune: Reached converged state.");
    if (ah_best(htask) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::tune: Error fetching best point.");
    }
    traversalOptionStr = ah_get_enum(htask, "traversalOption");
    traversalOption = TraversalOption::parseOptionExact(traversalOptionStr);
    dataLayoutOption = ah_get_enum(htask, "dataLayoutOption");
    newton3Option = ah_get_enum(htask, "newton3Option");
    cellSizeFactor = ah_get_real(htask, "cellSizeFactor");
    _currentConfig = Configuration(*compatibleTraversals::allCompatibleContainers(traversalOption).begin(),
                                   cellSizeFactor, traversalOption,
                                   DataLayoutOption::parseOptionExact(dataLayoutOption),
                                   Newton3Option::parseOptionExact(newton3Option));
    AutoPasLog(debug, "ActiveHarmony::tune: Selected optimal configuration {}.", _currentConfig.toString());
  }
  return !converged;
}

void ActiveHarmony::removeN3Option(Newton3Option option) {
  _allowedNewton3Options.erase(option);
  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
            "ActiveHarmony::removeN3Option: Removing all configurations with Newton 3 {} caused the search space to be empty!",
            option);
  }
}

const Configuration &ActiveHarmony::getCurrentConfiguration() const {
  return _currentConfig;
}

void ActiveHarmony::reset() {
  // free memory
  if (htask != nullptr) {
    ah_leave(htask);
    ah_kill(htask);
  }
  if (hdesc != nullptr) {
    ah_close(hdesc);
    ah_free(hdesc);
  }
  // initiate the tuning session
  hdesc = ah_alloc();
  if (hdesc == nullptr) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error allocating Harmony descriptor.");
  }

  if (ah_connect(hdesc, nullptr, 0) != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error connecting to Harmony session.");
  }
  // set up tuning parameters
  hdef_t *hdef = ah_def_alloc();
  if (hdef == nullptr) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error allocating definition descriptor");
  }

  if (ah_def_name(hdef, "AutoPas") != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error settings search name");
  }

  if (_allowedCellSizeFactors->isFinite() and _allowedCellSizeFactors->size() == 1) {} else {
    if (ah_def_real(hdef, "cellSizeFactor", _allowedCellSizeFactors->getMin(), _allowedCellSizeFactors->getMax(),
                    (_allowedCellSizeFactors->getMin(), _allowedCellSizeFactors->getMax()) / cellSizeSamples,
                    nullptr) !=
        0) {
      utils::ExceptionHandler::exception("ActiveHarmony::reset: Error defining real \"cellSizeFactor\"");
    }
  }

  if (ah_def_enum(hdef, "traversalOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error defining enum \"traversalOption\"");
  }
  for (auto &traversalOption : _allowedTraversalOptions) {
    if (ah_def_enum_value(hdef, "traversalOption", traversalOption.to_string().c_str()) != 0) {
      utils::ExceptionHandler::exception(
              "ActiveHarmony::reset: Error defining enum value for enum \"traversalOption\"");
    }
  }

  if (ah_def_enum(hdef, "dataLayoutOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error defining enum \"dataLayoutOption\"");
  }
  for (auto &dataLayoutOption : _allowedDataLayoutOptions) {
    if (ah_def_enum_value(hdef, "dataLayoutOption", dataLayoutOption.to_string().c_str()) != 0) {
      utils::ExceptionHandler::exception(
              "ActiveHarmony::reset: Error defining enum value for enum \"dataLayoutOption\"");
    }
  }

  if (ah_def_enum(hdef, "newton3Option", nullptr) != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error defining enum \"newton3Option\"");
  }
  for (auto &newton3Option : _allowedNewton3Options) {
    if (ah_def_enum_value(hdef, "newton3Option", newton3Option.to_string().c_str()) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::reset: Error defining enum value for enum \"newton3Option\"");
    }
  }

  ah_def_strategy(hdef, "nm.so");
  ah_def_cfg(hdef, "INIT_RADIUS", "0.5");
  // task initialization
  htask = ah_start(hdesc, hdef);
  ah_def_free(hdef);
  if (htask == nullptr) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error starting task.");
  }

  tune(false);
}

std::set<ContainerOption> ActiveHarmony::getAllowedContainerOptions() const {
  return _allowedContainerOptions;
}

bool ActiveHarmony::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _allowedContainerOptions.size() == 1 and
         (_allowedCellSizeFactors->isFinite() && _allowedCellSizeFactors->size() == 1) and
         _allowedTraversalOptions.size() == 1 and _allowedDataLayoutOptions.size() == 1 and
         _allowedNewton3Options.size() == 1;
}

bool ActiveHarmony::searchSpaceIsEmpty() const {
  return _allowedContainerOptions.empty() or
         (_allowedCellSizeFactors->isFinite() && _allowedCellSizeFactors->size() == 0) or
         _allowedTraversalOptions.empty() or _allowedDataLayoutOptions.empty() or _allowedNewton3Options.empty();
}

}  // namespace autopas
