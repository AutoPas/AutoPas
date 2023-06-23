/**
 * @file ActiveHarmony.cpp
 * @author Jakob Englhauser
 * @date 17.11.2022
 */

#include "ActiveHarmony.h"

autopas::ActiveHarmony::ActiveHarmony(const std::set<ContainerOption> &allowedContainerOptions,
                                      const autopas::NumberSet<double> &allowedCellSizeFactors,
                                      const std::set<TraversalOption> &allowedTraversalOptions,
                                      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                      const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                      const std::set<Newton3Option> &allowedNewton3Options,
                                      const autopas::MPIStrategyOption mpiStrategyOption,
                                      const autopas::AutoPas_MPI_Comm comm)
    : _allowedContainerOptions(allowedContainerOptions),
      _allowedCellSizeFactors(allowedCellSizeFactors.clone()),
      _allowedTraversalOptions(allowedTraversalOptions),
      _allowedLoadEstimatorOptions(allowedLoadEstimatorOptions),
      _allowedDataLayoutOptions(allowedDataLayoutOptions),
      _allowedNewton3Options(allowedNewton3Options),
      _currentConfig(),
      _mpiStrategyOption(mpiStrategyOption),
      _comm(comm),
      _nonLocalServer(getenv("HARMONY_HOST") != nullptr and mpiStrategyOption == MPIStrategyOption::divideAndConquer) {
  auto cellSizeDummy = NumberSetFinite<double>{-1};
  utils::AutoPasConfigurationCommunicator::distributeConfigurations(
      _allowedContainerOptions, cellSizeDummy, _allowedTraversalOptions, _allowedLoadEstimatorOptions,
      _allowedDataLayoutOptions, _allowedNewton3Options, 0, 1);

  AutoPasLog(DEBUG, "Possible container options: {}", autopas::utils::ArrayUtils::to_string(_allowedContainerOptions));
  AutoPasLog(DEBUG, "Possible traversal options: {}", autopas::utils::ArrayUtils::to_string(_allowedTraversalOptions));

  if (searchSpaceIsEmpty()) {
    autopas::utils::ExceptionHandler::exception("BayesianSearch: No valid configurations could be created.");
  }

  // set HARMONY_HOME environment variable; needed by active harmony library; the macro is set by cmake
  if (getenv("HARMONY_HOME") == nullptr) {
    putenv(const_cast<char *>(HARMONY_HOME));
  }

  reset(0);
}

autopas::ActiveHarmony::~ActiveHarmony() {
  if (htask != nullptr) {
    ah_leave(htask);
    ah_kill(htask);
  }
  if (hdesc != nullptr) {
    ah_close(hdesc);
    ah_free(hdesc);
  }
}

void autopas::ActiveHarmony::addEvidence(long time, size_t iteration) {
  if (searchSpaceIsTrivial() or searchSpaceIsEmpty()) {
    AutoPasLog(DEBUG, "ActiveHarmony::addEvidence: Search space is {}; did not report performance",
               searchSpaceIsTrivial() ? "trivial" : "empty");
  } else {
    auto perf = static_cast<double>(time);
    if (ah_report(htask, &perf) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::addEvidence: Error reporting performance to server");
    }
    _traversalTimes[_currentConfig] = time;
  }
}

template <class OptionClass>
OptionClass autopas::ActiveHarmony::fetchTuningParameter(const char *name, const std::set<OptionClass> &options) {
  OptionClass option;
  if (options.size() > 1) {
    option = decltype(option)::parseOptionExact(ah_get_enum(htask, name));
  } else if (options.size() == 1) {
    option = *options.begin();
  }
  return option;
}

void autopas::ActiveHarmony::fetchConfiguration() {
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

void autopas::ActiveHarmony::invalidateConfiguration() {
  auto worstPerf = std::numeric_limits<long>::max();
  addEvidence(worstPerf, 0);
  AutoPasLog(DEBUG, "ActiveHarmony::invalidateConfiguration: {}", _currentConfig.toString());
}

bool autopas::ActiveHarmony::tune(bool currentInvalid) {
  if (searchSpaceIsTrivial()) {
    fetchConfiguration();
    return false;
  } else if (searchSpaceIsEmpty()) {
    _currentConfig = Configuration();
    return false;
  }
  if (currentInvalid) {
    if (ah_converged(htask)) {
      AutoPasLog(DEBUG, "Active Harmony converged to invalid configuration; resetting active-harmony server.");
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
    if (_traversalTimes.find(_currentConfig) != _traversalTimes.end()) {
      // we already know the performance for this config
      addEvidence(static_cast<long>(_traversalTimes[_currentConfig]), 0);
      skipConfig = true;
    }

    auto converged = ah_converged(htask);
    if (converged) {
      // set configuration to optimum
      AutoPasLog(DEBUG, "ActiveHarmony::tune: Reached converged state.");
      if (ah_best(htask) != 0) {
        utils::ExceptionHandler::exception("ActiveHarmony::tune: Error fetching best point.");
      }
      fetchConfiguration();
      AutoPasLog(DEBUG, "ActiveHarmony::tune: Selected optimal configuration {}.", _currentConfig.toString());
      return false;
    }

    if (_nonLocalServer) {
      // When using a non-local server, it is possible that only tested configurations are fetched before the search
      // converges.
      // Because this is difficult to test for, the loop is simply ignored for non-local servers.
      return true;
    }
  } while (skipConfig);
  return true;
}

void autopas::ActiveHarmony::removeN3Option(Newton3Option option) {
  _allowedNewton3Options.erase(option);
  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "ActiveHarmony::removeN3Option: Removing all configurations with Newton 3 {} caused the search space to be "
        "empty!",
        option);
  }
  resetHarmony();
}

const autopas::Configuration &autopas::ActiveHarmony::getCurrentConfiguration() const { return _currentConfig; }

template <class OptionClass>
void autopas::ActiveHarmony::configureTuningParameter(hdef_t *hdef, const char *name,
                                                      const std::set<OptionClass> &options) {
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
    AutoPasLog(DEBUG, "ActiveHarmony::configureTuningParameter: Skipping trivial parameter {}", name);
  }
}

void autopas::ActiveHarmony::reset(size_t iteration) {
  _traversalTimes.clear();
  resetHarmony();
}

void autopas::ActiveHarmony::resetHarmony() {
  int rank, commSize;
  AutoPas_MPI_Comm_size(_comm, &commSize);
  AutoPas_MPI_Comm_rank(_comm, &rank);

  if (_mpiStrategyOption == MPIStrategyOption::divideAndConquer) {
    AutoPas_MPI_Barrier(_comm);
  }
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
    AutoPasLog(DEBUG, "Search space is {}; skipping harmony initialization.",
               searchSpaceIsTrivial() ? "trivial" : "empty");
    // Barrier in case the search space is bigger for some ranks
    if (_mpiStrategyOption == MPIStrategyOption::divideAndConquer) {
      AutoPas_MPI_Barrier(_comm);
    }
  } else {
    // initiate the tuning session
    hdesc = ah_alloc();
    if (hdesc == nullptr) {
      utils::ExceptionHandler::exception("ActiveHarmony::reset: Error allocating Harmony descriptor: {}", ah_error());
    }

    if (ah_connect(hdesc, nullptr, 0) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::reset: Error connecting to Harmony session: {}", ah_error());
    }

    if (not _nonLocalServer or rank == 0) {
      hdef_t *hdef = ah_def_alloc();
      if (hdef == nullptr) {
        utils::ExceptionHandler::exception("ActiveHarmony::reset: Error allocating definition descriptor: {}",
                                           ah_error());
      }
      setupTuningParameters(commSize, hdef);
      // task initialization
      htask = ah_start(hdesc, hdef);
      ah_def_free(hdef);
      if (htask == nullptr) {
        utils::ExceptionHandler::exception("ActiveHarmony::reset: Error starting task.");
      }

      if (_mpiStrategyOption == MPIStrategyOption::divideAndConquer) {
        AutoPas_MPI_Barrier(_comm);
      }
    } else {
      // only join a session if using the divideAndConquer mpi strategy, we are not rank 0 and a server is specified.
      if (_mpiStrategyOption == MPIStrategyOption::divideAndConquer) {
        AutoPas_MPI_Barrier(_comm);
      }

      // Everybody else may now join the master's new Harmony search.
      htask = ah_join(hdesc, "AutoPas");
      if (htask == nullptr) {
        utils::ExceptionHandler::exception("ActiveHarmony::reset: Error joining session");
      }
    }
  }

  tune(false);
}

std::set<autopas::ContainerOption> autopas::ActiveHarmony::getAllowedContainerOptions() const {
  return _allowedContainerOptions;
}

bool autopas::ActiveHarmony::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _allowedContainerOptions.size() == 1 and
         (_allowedCellSizeFactors->isFinite() and _allowedCellSizeFactors->size() == 1) and
         _allowedTraversalOptions.size() == 1 and _allowedDataLayoutOptions.size() == 1 and
         _allowedNewton3Options.size() == 1;
}

bool autopas::ActiveHarmony::searchSpaceIsEmpty() const {
  return _allowedContainerOptions.empty() or
         (_allowedCellSizeFactors->isFinite() and _allowedCellSizeFactors->size() == 0) or
         _allowedTraversalOptions.empty() or _allowedDataLayoutOptions.empty() or _allowedNewton3Options.empty();
}

void autopas::ActiveHarmony::setupTuningParameters(int commSize, hdef_t *hdef) {
  if (ah_def_name(hdef, "AutoPas") != 0) {
    utils::ExceptionHandler::exception("ActiveHarmony::reset: Error setting search name: {}", ah_error());
  }

  if (_allowedCellSizeFactors->isFinite()) {  // finite cell-size factors => define parameter as enum
    if (_allowedCellSizeFactors->size() == 1) {
      AutoPasLog(DEBUG, "ActiveHarmony::reset: Skipping trivial parameter {}", cellSizeFactorsName);
    } else if (_allowedCellSizeFactors->size() > 1) {
      AutoPasLog(DEBUG, "ActiveHarmony::reset: Finite cell-size factors; defining parameter as enum");
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
    AutoPasLog(DEBUG, "ActiveHarmony::reset: Infinite cell-size factors; defining parameter as real");
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
  configureTuningParameter(hdef, loadEstimatorOptionName, _allowedLoadEstimatorOptions);

  // use ActiveHarmony's implementation of the Nelder-Mead method
  ah_def_strategy(hdef, "pro.so");
  // set the size of the initial simplex (as portion of the total search space)
  ah_def_cfg(hdef, "INIT_RADIUS", "0.7");

  if (_mpiStrategyOption == MPIStrategyOption::divideAndConquer) {
    // set the size of the tuning session
    char numbuf[12];
    snprintf(numbuf, sizeof(numbuf), "%d", commSize);
    ah_def_cfg(hdef, "CLIENT_COUNT", numbuf);
  }
}
long autopas::ActiveHarmony::getEvidence(autopas::Configuration configuration) const {
  return static_cast<long>(_traversalTimes.at(configuration));
}

bool autopas::ActiveHarmony::smoothedHomogeneityAndMaxDensityNeeded() const { return false; }
