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
   * Constructor. Note that ActiveHarmony assumes every traversal option is only applicable for one container.
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
    // set HARMONY_HOME environment variable; needed by active harmony library; the macro is set by cmake
    if (getenv("HARMONY_HOME") == nullptr) {
      putenv(const_cast<char *>(HARMONY_HOME));
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
  /**
   * Pointer for the connection to the ActiveHarmony server.
   */
  hdesc_t *hdesc;
  /**
   * Pointer to the ActiveHarmony tuning task defining the tuning parameters and tuning process.
   */
  htask_t *htask;

  std::set<ContainerOption> _allowedContainerOptions;
  std::unique_ptr<NumberSet < double>> _allowedCellSizeFactors;
  std::set<TraversalOption> _allowedTraversalOptions;
  std::set<DataLayoutOption> _allowedDataLayoutOptions;
  std::set<Newton3Option> _allowedNewton3Options;

  Configuration _currentConfig;

  /**
   * Fetch parameter-values from harmony server and update _currentConfig.
   */
  inline void fetchConfiguration();

  /**
   * Invalidates the current configuration by reporting the worst possible performance to the harmony server.
   * Always returns true to allow being used in boolean statements.
   * @return true
   */
  inline bool invalidateConfiguration();

  /**
   * Define the tuning parameter in the ActiveHarmony tuning definition as enum and set possible values.
   * @tparam O Option class.
   * @param hdef Pointer to ActiveHarmony tuning definition.
   * @param name Name of the tuning parameter.
   * @param options Set of possible values of the tuning parameter.
   */
  template <class O>
  inline void configureTuningParameter(hdef_t *hdef, const char *name, std::set<O> options);

  static constexpr int cellSizeSamples = 100;

  static constexpr const char *traversalOptionName = "traversalOption";
  static constexpr const char *dataLayoutOptionName = "dataLayoutOption";
  static constexpr const char *cellSizeFactorName = "cellSizeFactor";
  static constexpr const char *newton3OptionName = "newton3Option";

};

void ActiveHarmony::addEvidence(long time) {
  auto perf = (double) time;
  if (ah_report(htask, &perf) != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::addEvidence: Error reporting performance to server");
  }
}

void ActiveHarmony::fetchConfiguration() {
  TraversalOption traversalOption;
  if (_allowedTraversalOptions.size() > 1) {
    traversalOption = TraversalOption::parseOptionExact(ah_get_enum(htask, traversalOptionName));
  } else if (_allowedTraversalOptions.size() == 1) {
    traversalOption = *_allowedTraversalOptions.begin();
  }

  DataLayoutOption dataLayoutOption;
  if (_allowedDataLayoutOptions.size() > 1) {
    dataLayoutOption = DataLayoutOption::parseOptionExact(ah_get_enum(htask, dataLayoutOptionName));
  } else if (_allowedDataLayoutOptions.size() == 1) {
    dataLayoutOption = *_allowedDataLayoutOptions.begin();
  }

  Newton3Option newton3Option;
  if (_allowedNewton3Options.size() > 1) {
    newton3Option = Newton3Option::parseOptionExact(ah_get_enum(htask, newton3OptionName));
  } else if (_allowedNewton3Options.size() == 1) {
    newton3Option = *_allowedNewton3Options.begin();
  }

  double cellSizeFactor;
  if (_allowedCellSizeFactors->isFinite() and _allowedCellSizeFactors->size() == 1) {
    cellSizeFactor = _allowedCellSizeFactors->getMin();
  } else {
    cellSizeFactor = ah_get_real(htask, cellSizeFactorName);
  }

  _currentConfig = Configuration(*compatibleTraversals::allCompatibleContainers(traversalOption).begin(),
                                 cellSizeFactor, traversalOption, dataLayoutOption, newton3Option);
}

bool ActiveHarmony::invalidateConfiguration() {
  auto worstPerf = std::numeric_limits<double>::max();
  addEvidence(worstPerf);
  return true;
}

bool ActiveHarmony::tune(bool currentInvalid) {
  if (currentInvalid) {
    invalidateConfiguration();
  }

  // get configurations from server until valid newton3 option is found
  do {
    if (ah_fetch(htask) < 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::tune: Error fetching values from server");
    }
    fetchConfiguration();
  } while (_allowedNewton3Options.find(_currentConfig.newton3) == _allowedNewton3Options.end() and
           invalidateConfiguration()); // short-circuit evaluation makes this only execute if newton3 option is invalid

  auto converged = ah_converged(htask);
  if (converged) {
    // set configuration to optimum
    AutoPasLog(debug, "ActiveHarmony::tune: Reached converged state.");
    if (ah_best(htask) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::tune: Error fetching best point.");
    }
    fetchConfiguration();
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

template<class O>
void ActiveHarmony::configureTuningParameter(hdef_t *hdef, const char *name, std::set<O> options) {
  if (options.size() > 1) {
    if (ah_def_enum(hdef, name, nullptr) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::configureTuningParameter: Error defining enum \"{}\"", name);
    }
    for (auto &option : options) {
      if (ah_def_enum_value(hdef, name, option.to_string().c_str()) != 0) {
        utils::ExceptionHandler::exception(
                "ActiveHarmony::configureTuningParameter: Error defining enum value for enum \"{}\"", name);
      }
    }
  } else {
    AutoPasLog(debug, "ActiveHarmony::configureTuningParameter: Skipping trivial parameter {}", name);
  }
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

  if (not _allowedCellSizeFactors->isFinite() or _allowedCellSizeFactors->size() > 1) {
    if (ah_def_real(hdef, cellSizeFactorName, _allowedCellSizeFactors->getMin(), _allowedCellSizeFactors->getMax(),
                    (_allowedCellSizeFactors->getMin(), _allowedCellSizeFactors->getMax()) / cellSizeSamples,
                    nullptr) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::reset: Error defining real \"{}\"", cellSizeFactorName);
    }
  } else {
    AutoPasLog(debug, "ActiveHarmony::reset: Skipping trivial parameter {}", cellSizeFactorName);
  }

  configureTuningParameter(hdef, traversalOptionName, _allowedTraversalOptions);
  configureTuningParameter(hdef, dataLayoutOptionName, _allowedDataLayoutOptions);
  configureTuningParameter(hdef, newton3OptionName, _allowedNewton3Options);

  // use ActiveHarmony's implementation of the Nelder-Mead method
  ah_def_strategy(hdef, "nm.so");
  // set the size of the initial simplex (as portion of the total search space)
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
         (_allowedCellSizeFactors->isFinite() and _allowedCellSizeFactors->size() == 1) and
         _allowedTraversalOptions.size() == 1 and _allowedDataLayoutOptions.size() == 1 and
         _allowedNewton3Options.size() == 1;
}

bool ActiveHarmony::searchSpaceIsEmpty() const {
  return _allowedContainerOptions.empty() or
         (_allowedCellSizeFactors->isFinite() and _allowedCellSizeFactors->size() == 0) or
         _allowedTraversalOptions.empty() or _allowedDataLayoutOptions.empty() or _allowedNewton3Options.empty();
}

}  // namespace autopas
