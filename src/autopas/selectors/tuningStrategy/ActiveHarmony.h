/**
 * @file ActiveHarmony.h
 * @author Jakob Englhauser
 * @date 21.10.19
 */

#pragma once

#include <string>
#include <sstream>
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
  ActiveHarmony(const std::set<ContainerOption> &allowedContainerOptions,
                const NumberSetFinite<double> &allowedCellSizeFactors,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options)
      : _allowedContainerOptions(allowedContainerOptions),
        _allowedCellSizeFactors(allowedCellSizeFactors.clone()),
        _allowedTraversalOptions(allowedTraversalOptions),
        _allowedDataLayoutOptions(allDataLayoutOptions),
        _allowedNewton3Options(allNewton3Options) {

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

  bool searchSpaceIsTrivial() const override;

  bool searchSpaceIsEmpty() const override;

  const Configuration &getCurrentConfiguration() const override;

  void reset() override;

  std::set<ContainerOption> getAllowedContainerOptions() const override;

 private:
  hdesc_t *hdesc;
  htask_t *htask;

  std::set<ContainerOption> _allowedContainerOptions;
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors; // Maybe use min, max, stepsize
  std::set<TraversalOption> _allowedTraversalOptions;
  std::set<DataLayoutOption> _allowedDataLayoutOptions;
  std::set<Newton3Option> _allowedNewton3Options;
};

void ActiveHarmony::addEvidence(long time) {
  auto perf = (double) time;
  if (ah_report(htask, &perf) != 0) {
    utils::ExceptionHandler::exception("Error reporting performance to server");
  }
}

bool ActiveHarmony::tune(bool currentInvalid) {
  if (ah_fetch(htask) < 0) {
    utils::ExceptionHandler::exception("Error fetching values from server");
  }
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
/*
const Configuration &ActiveHarmony::getCurrentConfiguration() const {
  // TODO map tuning parameters to config
  return ;
}*/

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
  // TODO check for valid container, traversal combination
  /*
  if (ah_def_enum(hdef, "containerOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"containerOption\"");
  }
  for (auto &containerOption : _allowedContainerOptions) {
    std::stringstream ss;
    ss<<containerOption; // convert enum to string
    if (ah_def_enum_value(hdef, "containerOption", ss.str().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"containerOption\"");
    }
  }

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

  if (ah_def_enum(hdef, "traversalOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"traversalOption\"");
  }
  for (auto &traversalOption : _allowedTraversalOptions) {
    std::stringstream ss;
    ss<<traversalOption; // convert enum to string
    if (ah_def_enum_value(hdef, "traversalOption", ss.str().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"traversalOption\"");
    }
  }

  if (ah_def_enum(hdef, "dataLayoutOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"dataLayoutOption\"");
  }
  for (auto &dataLayoutOption : _allowedDataLayoutOptions) {
    std::stringstream ss;
    ss<<dataLayoutOption; // convert enum to string
    if (ah_def_enum_value(hdef, "dataLayoutOption", ss.str().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"dataLayoutOption\"");
    }
  }

  if (ah_def_enum(hdef, "containerOption", nullptr) != 0) {
    utils::ExceptionHandler::exception("Error defining enum \"containerOption\"");
  }
  for (auto &newton3Option : _allowedNewton3Options) {
    std::stringstream ss;
    ss<<newton3Option; // convert enum to string
    if (ah_def_enum_value(hdef, "newton3Option", ss.str().c_str()) != 0) {
      utils::ExceptionHandler::exception("Error defining enum value for enum \"newton3Option\"");
    }
  }*/

  // task configuration TODO
  ah_def_strategy(hdef, "pro.so");
  ah_def_cfg(hdef, "INIT_RADIUS", "0.5");
  ah_def_layers(hdef, "log.so:agg.so");
  ah_def_cfg(hdef, "LOG_FILE", "/tmp/tuning.run");
  ah_def_cfg(hdef, "AGG_TIMES", "10");
  ah_def_cfg(hdef, "AGG_FUNC", "median");
  // task initialization
  // TODO bind tuning parameters to local variables, or maybe retrieve parameters in getConfiguration ?
  htask = ah_start(hdesc, hdef);
  ah_def_free(hdef);
  htask = ah_start(hdesc, hdef);
}

std::set<ContainerOption> ActiveHarmony::getAllowedContainerOptions() const {
  return _allowedContainerOptions;
}

bool ActiveHarmony::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _allowedContainerOptions.size() == 1 and (_allowedCellSizeFactors->isFinite() && _allowedCellSizeFactors->size() == 1) and
         _allowedTraversalOptions.size() == 1 and _allowedDataLayoutOptions.size() == 1 and _allowedNewton3Options.size() == 1;
}

bool ActiveHarmony::searchSpaceIsEmpty() const {
  return _allowedContainerOptions.empty() or (_allowedCellSizeFactors->isFinite() && _allowedCellSizeFactors->size() == 0) or
         _allowedTraversalOptions.empty() or _allowedDataLayoutOptions.empty() or _allowedNewton3Options.empty();
}

}  // namespace autopas
