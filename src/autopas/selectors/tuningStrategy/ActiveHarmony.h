/**
 * @file ActiveHarmony.h
 * @author Jakob Englhauser
 * @date 21.10.19
 */

#pragma once

#include <string>
#include <autopas/selectors/FeatureVector.h>
#include <autopas/containers/CompatibleTraversals.h>
#include <cfloat>
#include "TuningStrategyInterface.h"
#include "hclient.h"

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
  ActiveHarmony(const NumberSetFinite<double> &allowedCellSizeFactors,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options)
      : _allowedCellSizeFactors(allowedCellSizeFactors.clone()),
        _allowedTraversalOptions(allowedTraversalOptions),
        _allowedDataLayoutOptions(allowedDataLayoutOptions),
        _allowedNewton3Options(allowedNewton3Options),
        _currentConfig() {

    reset();
  }

  ~ActiveHarmony() override {
    if (htask != nullptr) {
      ah_leave(htask);
    }
    if (hdesc != nullptr) {
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

  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors; // Maybe use min, max, stepsize
  std::set<TraversalOption> _allowedTraversalOptions;
  std::set<DataLayoutOption> _allowedDataLayoutOptions;
  std::set<Newton3Option> _allowedNewton3Options;

  FeatureVector _currentConfig;
};

void ActiveHarmony::addEvidence(long time) {
  auto perf = (double) time;
  if (ah_report(htask, &perf) != 0) {
    utils::ExceptionHandler::exception("Error reporting performance to server");
  }
}

bool ActiveHarmony::tune(bool currentInvalid) {
  if (currentInvalid) {
    auto perf = DBL_MAX;
    ah_report(htask, &perf);
  }
  if (ah_fetch(htask) < 0) {
    utils::ExceptionHandler::exception("Error fetching values from server");
  }
  std::string traversalOptionStr = ah_get_enum(htask, "traversalOption");
  std::transform(traversalOptionStr.begin(), traversalOptionStr.end(), traversalOptionStr.begin(), ::tolower);
  auto traversalOption = *TraversalOption::parseOptions(traversalOptionStr).begin();
  std::string dataLayoutOption = ah_get_enum(htask, "dataLayoutOption");
  std::transform(dataLayoutOption.begin(), dataLayoutOption.end(), dataLayoutOption.begin(), ::tolower);
  std::string newton3Option = ah_get_enum(htask, "newton3Option");
  std::transform(newton3Option.begin(), newton3Option.end(), newton3Option.begin(), ::tolower);
  _currentConfig = FeatureVector(*compatibleTraversals::allCompatibleContainers(traversalOption).begin(), 1.0, traversalOption, *DataLayoutOption::parseOptions(dataLayoutOption).begin(), *Newton3Option::parseOptions(newton3Option).begin());
  return !ah_converged(htask);
}

void ActiveHarmony::removeN3Option(Newton3Option option) {
  _allowedNewton3Options.erase(option);
  reset();
  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
            "Removing all configurations with Newton 3 {} caused the search space to be empty!", option);
  }
}

const Configuration &ActiveHarmony::getCurrentConfiguration() const {
  return _currentConfig;
}

void ActiveHarmony::reset() {
// initiate the tuning session
  hdesc = ah_alloc();
  if (hdesc == nullptr) {
    utils::ExceptionHandler::exception("Error allocating Harmony descriptor.");
    // TODO handle error
  }

  if (ah_connect(hdesc, nullptr, 0) != 0) {
    utils::ExceptionHandler::exception("Error connecting to Harmony session.");
  }

  hdef_t *hdef = ah_def_alloc();
  if (hdef == nullptr) {
    utils::ExceptionHandler::exception("Error allocating definition descriptor");
  }

  if (ah_def_name(hdef, "AutoPas") != 0) {
    utils::ExceptionHandler::exception("Error settings search name");
  }
  // tuning parameters


/*
  if (ah_def_enum(hdef, "cellSizeFactor", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"cellSizeFactor\"");
  }
  for (auto &cellSizeFactor : _allowedCellSizeFactors->getAll()) {
    std::stringstream ss;
    ss<<cellSizeFactor; // convert enum to string
    if (ah_def_enum_value(hdef, "cellSizeFactor", ss.str().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"cellSizeFactor\"");
    }
  }
*/
  if (ah_def_enum(hdef, "traversalOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"traversalOption\"");
  }
  for (auto &traversalOption : _allowedTraversalOptions) {
    if (ah_def_enum_value(hdef, "traversalOption", traversalOption.to_string().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"traversalOption\"");
    }
  }

  if (ah_def_enum(hdef, "dataLayoutOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"dataLayoutOption\"");
  }
  for (auto &dataLayoutOption : _allowedDataLayoutOptions) {
    if (ah_def_enum_value(hdef, "dataLayoutOption", dataLayoutOption.to_string().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"dataLayoutOption\"");
    }
  }

  if (ah_def_enum(hdef, "newton3Option", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"newton3Option\"");
  }
  for (auto &newton3Option : _allowedNewton3Options) {
    if (ah_def_enum_value(hdef, "newton3Option", newton3Option.to_string().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"newton3Option\"");
    }
  }

  // task configuration TODO
  ah_def_strategy(hdef, "pro.so");
  ah_def_cfg(hdef, "INIT_RADIUS", "0.5");
  ah_def_layers(hdef, "log.so");
  ah_def_cfg(hdef, "LOG_FILE", "/tmp/tuning.run");
  // task initialization
  // TODO bind tuning parameters to local variables, or maybe retrieve parameters in getConfiguration ?
  htask = ah_start(hdesc, hdef);
  ah_def_free(hdef);

  _currentConfig = FeatureVector(*compatibleTraversals::allCompatibleContainers(*_allowedTraversalOptions.begin()).begin(), 1.0, *_allowedTraversalOptions.begin(), *_allowedDataLayoutOptions.begin(), *_allowedNewton3Options.begin());
}

std::set<ContainerOption> ActiveHarmony::getAllowedContainerOptions() const {
  std::set<ContainerOption> allowedContainerOptions;
  for (auto traversalOption : _allowedTraversalOptions) {
    auto compatibleContainerOptions = compatibleTraversals::allCompatibleContainers(traversalOption);
    allowedContainerOptions.insert(compatibleContainerOptions.begin(), compatibleContainerOptions.end());
  }
  return allowedContainerOptions;
}

bool ActiveHarmony::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return (_allowedCellSizeFactors->isFinite() && _allowedCellSizeFactors->size() == 1) and
         _allowedTraversalOptions.size() == 1 and _allowedDataLayoutOptions.size() == 1 and _allowedNewton3Options.size() == 1;
}

bool ActiveHarmony::searchSpaceIsEmpty() const {
  return (_allowedCellSizeFactors->isFinite() && _allowedCellSizeFactors->size() == 0) or
         _allowedTraversalOptions.empty() or _allowedDataLayoutOptions.empty() or _allowedNewton3Options.empty();
}

}  // namespace autopas
