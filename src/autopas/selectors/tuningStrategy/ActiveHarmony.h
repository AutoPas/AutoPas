/**
 * @file ActiveHarmony.h
 * @author Jakob Englhauser
 * @date 21.10.2019
 */

#pragma once

#include <string>

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/loadEstimators.h"
#include "autopas/selectors/Configuration.h"
#include "hclient.h"

namespace autopas {

/**
 * Interface to the Active Harmony tuning framework.
 */
class ActiveHarmony : public TuningStrategyInterface {
 public:
  /**
   * Constructor. Note that ActiveHarmony assumes every traversal option is only applicable for one container.
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  ActiveHarmony(const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
                const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
                const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
                const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
                const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
                const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions())
      : _allowedContainerOptions(),
        _allowedCellSizeFactors(allowedCellSizeFactors.clone()),
        _allowedTraversalOptions(allowedTraversalOptions),
        _allowedLoadEstimatorOptions(allowedLoadEstimatorOptions),
        _allowedDataLayoutOptions(allowedDataLayoutOptions),
        _allowedNewton3Options(allowedNewton3Options),
        _currentConfig() {
    AutoPasLog(debug, "Reducing traversal and container options to possible combinations");
    for (auto &traversalOption : _allowedTraversalOptions) {
      auto compatibleContainers = compatibleTraversals::allCompatibleContainers(traversalOption);
      if (compatibleContainers.empty() or
          allowedContainerOptions.find(*compatibleContainers.begin()) == allowedContainerOptions.end()) {
        _allowedTraversalOptions.erase(traversalOption);
      } else {
        _allowedContainerOptions.emplace(*compatibleContainers.begin());
      }
    }
    AutoPasLog(debug, "Possible container options: {}",
               autopas::utils::ArrayUtils::to_string(_allowedContainerOptions));
    AutoPasLog(debug, "Possible traversal options: {}",
               autopas::utils::ArrayUtils::to_string(_allowedTraversalOptions));

    // set HARMONY_HOME environment variable; needed by active harmony library; the macro is set by cmake
    if (getenv("HARMONY_HOME") == nullptr) {
      putenv(const_cast<char *>(HARMONY_HOME));
    }

    reset(0);
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

  inline void addEvidence(long time, size_t iteration) override;

  inline bool tune(bool currentInvalid) override;

  inline void removeN3Option(Newton3Option option) override;

  inline bool searchSpaceIsTrivial() const override;

  inline bool searchSpaceIsEmpty() const override;

  inline const Configuration &getCurrentConfiguration() const override;

  inline void reset(size_t iteration) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override;

 private:
  /**
   * Pointer for the connection to the ActiveHarmony server.
   */
  hdesc_t *hdesc = nullptr;
  /**
   * Pointer to the ActiveHarmony tuning task defining the tuning parameters and tuning process.
   */
  htask_t *htask = nullptr;

  std::set<ContainerOption> _allowedContainerOptions;
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors;
  std::set<TraversalOption> _allowedTraversalOptions;
  std::set<LoadEstimatorOption> _allowedLoadEstimatorOptions;
  std::set<DataLayoutOption> _allowedDataLayoutOptions;
  std::set<Newton3Option> _allowedNewton3Options;

  Configuration _currentConfig;

  /**
   * Traversal times for configurations. Used to save evidence when resetting active-harmony server.
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;

  /**
   * Resets the Harmony Server but keeps evidence.
   */
  inline void resetHarmony();

  /**
   * Fetch parameter-values from harmony server and update _currentConfig.
   */
  inline void fetchConfiguration();

  /**
   * Invalidates the current configuration by reporting the worst possible performance to the harmony server.
   */
  inline void invalidateConfiguration();

  /**
   * Define the tuning parameter in the ActiveHarmony tuning definition as enum and set possible values.
   * @tparam O Option class.
   * @param hdef Pointer to ActiveHarmony tuning definition.
   * @param name Name of the tuning parameter.
   * @param options Set of possible values of the tuning parameter.
   */
  template <class O>
  inline void configureTuningParameter(hdef_t *hdef, const char *name, const std::set<O> options);

  /**
   * Fetch value for enum-type tuning parameter.
   * @tparam O Option class.
   * @param name Name of the tuning parameter.
   * @param options Set of all allowed options.
   * @return Value for tuning parameter.
   */
  template <class O>
  inline O fetchTuningParameter(const char *name, const std::set<O> options);

  static constexpr int cellSizeSamples = 100;

  static constexpr const char *traversalOptionName = "traversalOption";
  static constexpr const char *loadEstimatorOptionName = "loadEstimatorOption";
  static constexpr const char *dataLayoutOptionName = "dataLayoutOption";
  static constexpr const char *cellSizeFactorsName = "cellSizeFactor";
  static constexpr const char *newton3OptionName = "newton3Option";
};

void ActiveHarmony::addEvidence(long time, size_t iteration) {
  if (searchSpaceIsTrivial() or searchSpaceIsEmpty()) {
    AutoPasLog(debug, "ActiveHarmony::addEvidence: Search space is {}; did not report performance",
               searchSpaceIsTrivial() ? "trivial" : "empty");
  } else {
    auto perf = (double)time;
    if (ah_report(htask, &perf) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::addEvidence: Error reporting performance to server");
    }
    _traversalTimes[_currentConfig] = time;
  }
}

template <class O>
O ActiveHarmony::fetchTuningParameter(const char *name, const std::set<O> options) {
  O option;
  if (options.size() > 1) {
    option = decltype(option)::parseOptionExact(ah_get_enum(htask, name));
  } else if (options.size() == 1) {
    option = *options.begin();
  }
  return option;
}

void ActiveHarmony::fetchConfiguration() {
  TraversalOption traversalOption = fetchTuningParameter(traversalOptionName, _allowedTraversalOptions);
  DataLayoutOption dataLayoutOption = fetchTuningParameter(dataLayoutOptionName, _allowedDataLayoutOptions);
  Newton3Option newton3Option = fetchTuningParameter(newton3OptionName, _allowedNewton3Options);
  ContainerOption containerOption = *compatibleTraversals::allCompatibleContainers(traversalOption).begin();
  auto applicableLoadEstimators =
      loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, _allowedLoadEstimatorOptions);
  LoadEstimatorOption loadEstimatorOption = fetchTuningParameter(loadEstimatorOptionName, applicableLoadEstimators);

  double cellSizeFactor = 0;
  if (_allowedCellSizeFactors->isFinite()) {
    if (_allowedCellSizeFactors->size() == 1) {
      cellSizeFactor = _allowedCellSizeFactors->getMin();
    } else if (_allowedCellSizeFactors->size() > 1) {
      cellSizeFactor = std::strtod(ah_get_enum(htask, cellSizeFactorsName), nullptr);
    }
  } else {
    cellSizeFactor = ah_get_real(htask, cellSizeFactorsName);
  }

  _currentConfig = Configuration(containerOption, cellSizeFactor, traversalOption, loadEstimatorOption,
                                 dataLayoutOption, newton3Option);
}

void ActiveHarmony::invalidateConfiguration() {
  auto worstPerf = std::numeric_limits<double>::max();
  addEvidence(worstPerf, 0);
}

bool ActiveHarmony::tune(bool currentInvalid) {
  if (searchSpaceIsTrivial()) {
    fetchConfiguration();
    return false;
  } else if (searchSpaceIsEmpty()) {
    _currentConfig = Configuration();
    return false;
  }
  if (currentInvalid) {
    if (ah_converged(htask)) {
      AutoPasLog(debug, "Active Harmony converged to invalid configuration; resetting active-harmony server.");
      resetHarmony();
    } else {
      invalidateConfiguration();
    }
  }

  // get configurations from server until new configuration with valid newton3 option is found
  bool skipConfig;
  do {
    skipConfig = false;
    if (ah_fetch(htask) < 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::tune: Error fetching values from server");
    }
    fetchConfiguration();
    if (_traversalTimes.find(_currentConfig) != _traversalTimes.end()) {  // we already know performance for this config
      addEvidence(_traversalTimes[_currentConfig], 0);
      skipConfig = true;
    }

    auto converged = ah_converged(htask);
    if (converged) {
      // set configuration to optimum
      AutoPasLog(debug, "ActiveHarmony::tune: Reached converged state.");
      if (ah_best(htask) != 0) {
        utils::ExceptionHandler::exception("ActiveHarmony::tune: Error fetching best point.");
      }
      fetchConfiguration();
      AutoPasLog(debug, "ActiveHarmony::tune: Selected optimal configuration {}.", _currentConfig.toString());
      return false;
    }
  } while (skipConfig);
  return true;
}

void ActiveHarmony::removeN3Option(Newton3Option option) {
  _allowedNewton3Options.erase(option);
  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "ActiveHarmony::removeN3Option: Removing all configurations with Newton 3 {} caused the search space to be "
        "empty!",
        option);
  }
  resetHarmony();
}

const Configuration &ActiveHarmony::getCurrentConfiguration() const { return _currentConfig; }

template <class O>
void ActiveHarmony::configureTuningParameter(hdef_t *hdef, const char *name, const std::set<O> options) {
  if (options.size() > 1) {                       // only parameters with more than 1 possible options should be tuned
    if (ah_def_enum(hdef, name, nullptr) != 0) {  // define parameter
      utils::ExceptionHandler::exception("ActiveHarmony::configureTuningParameter: Error defining enum \"{}\"", name);
    }
    for (auto &option : options) {  // define possible values for parameter
      if (ah_def_enum_value(hdef, name, option.to_string().c_str()) != 0) {
        utils::ExceptionHandler::exception(
            "ActiveHarmony::configureTuningParameter: Error defining enum value for enum \"{}\"", name);
      }
    }
  } else {
    AutoPasLog(debug, "ActiveHarmony::configureTuningParameter: Skipping trivial parameter {}", name);
  }
}

void ActiveHarmony::reset(size_t iteration) {
  _traversalTimes.clear();
  resetHarmony();
}

void ActiveHarmony::resetHarmony() {
  // free memory
  if (htask != nullptr) {
    ah_leave(htask);
    ah_kill(htask);
  }
  if (hdesc != nullptr) {
    ah_close(hdesc);
    ah_free(hdesc);
  }

  if (searchSpaceIsTrivial() or searchSpaceIsEmpty()) {
    AutoPasLog(debug, "Search space is {}; skipping harmony initialization.",
               searchSpaceIsTrivial() ? "trivial" : "empty");
  } else {
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
      utils::ExceptionHandler::exception("ActiveHarmony::reset: Error setting search name");
    }

    if (_allowedCellSizeFactors->isFinite()) {  // finite cell-size factors => define parameter as enum
      if (_allowedCellSizeFactors->size() == 1) {
        AutoPasLog(debug, "ActiveHarmony::reset: Skipping trivial parameter {}", cellSizeFactorsName);
      } else if (_allowedCellSizeFactors->size() > 1) {
        AutoPasLog(debug, "ActiveHarmony::reset: Finite cell-size factors; defining parameter as enum");
        if (ah_def_enum(hdef, cellSizeFactorsName, nullptr) != 0) {
          utils::ExceptionHandler::exception("ActiveHarmony::configureTuningParameter: Error defining enum \"{}\"",
                                             cellSizeFactorsName);
        }
        for (auto cellSizeFactor : _allowedCellSizeFactors->getAll()) {
          if (ah_def_enum_value(hdef, cellSizeFactorsName, std::to_string(cellSizeFactor).c_str()) != 0) {
            utils::ExceptionHandler::exception(
                "ActiveHarmony::configureTuningParameter: Error defining enum value for enum \"{}\"",
                cellSizeFactorsName);
          }
        }
      }
    } else {  // infinite cell-size factors => define parameter as real
      AutoPasLog(debug, "ActiveHarmony::reset: Infinite cell-size factors; defining parameter as real");
      if (ah_def_real(hdef, cellSizeFactorsName, _allowedCellSizeFactors->getMin(), _allowedCellSizeFactors->getMax(),
                      (_allowedCellSizeFactors->getMax() - _allowedCellSizeFactors->getMin()) / cellSizeSamples,
                      nullptr) != 0) {
        utils::ExceptionHandler::exception("ActiveHarmony::reset: Error defining real \"{}\"", cellSizeFactorsName);
      }
    }
    // set up other parameters
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
  }

  tune(false);
}

std::set<ContainerOption> ActiveHarmony::getAllowedContainerOptions() const { return _allowedContainerOptions; }

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
