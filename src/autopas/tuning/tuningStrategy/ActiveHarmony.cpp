/**
 * @file ActiveHarmony.cpp
 * @author Jakob Englhauser
 * @date 17.11.2022
 */

#include "ActiveHarmony.h"

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"

namespace autopas {

ActiveHarmony::ActiveHarmony(const InteractionTypeOption &interactionType,
                             const std::set<ContainerOption> &allowedContainerOptions,
                             const NumberSet<double> &allowedCellSizeFactors,
                             const std::set<TraversalOption> &allowedTraversalOptions,
                             const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                             const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                             const std::set<Newton3Option> &allowedNewton3Options, bool mpiDivideAndConquer,
                             AutoPas_MPI_Comm comm)
    : _interactionType(interactionType),
      _allowedContainerOptions(allowedContainerOptions),
      _allowedCellSizeFactors(allowedCellSizeFactors.clone()),
      _allowedTraversalOptions(allowedTraversalOptions),
      _allowedLoadEstimatorOptions(allowedLoadEstimatorOptions),
      _allowedDataLayoutOptions(allowedDataLayoutOptions),
      _allowedNewton3Options(allowedNewton3Options),
      _mpiDivideAndConquer(mpiDivideAndConquer),
      _comm(comm),
      _nonLocalServer(getenv("HARMONY_HOST") != nullptr and mpiDivideAndConquer) {
#ifdef AUTOPAS_ENABLE_HARMONY
  if (searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception("ActiveHarmony: No valid configurations could be created.");
  }

  // set HARMONY_HOME environment variable; needed by active harmony library; the macro is set by cmake
  if (getenv("HARMONY_HOME") == nullptr) {
    putenv(const_cast<char *>(HARMONY_HOME));
  }
#endif
}

ActiveHarmony::~ActiveHarmony() {
#ifdef AUTOPAS_ENABLE_HARMONY
  if (htask != nullptr) {
    ah_leave(htask);
    ah_kill(htask);
  }
  if (hdesc != nullptr) {
    ah_close(hdesc);
    ah_free(hdesc);
  }
#endif
}

void ActiveHarmony::addEvidence(const Configuration &configuration, const Evidence &evidence) {
#ifdef AUTOPAS_ENABLE_HARMONY
  if (searchSpaceIsTrivial() or searchSpaceIsEmpty()) {
    AutoPasLog(DEBUG, "ActiveHarmony::addEvidence: Search space is {}; did not report performance",
               searchSpaceIsTrivial() ? "trivial" : "empty");
  } else {
    auto perf = static_cast<double>(evidence.value);
    if (ah_report(htask, &perf) != 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::addEvidence: Error reporting performance to server");
    }
  }
#endif
}

template <class OptionClass>
OptionClass ActiveHarmony::fetchTuningParameter(const char *name, const std::set<OptionClass> &options) {
#ifdef AUTOPAS_ENABLE_HARMONY
  OptionClass option;
  if (options.size() > 1) {
    option = decltype(option)::parseOptionExact(ah_get_enum(htask, name));
  } else if (options.size() == 1) {
    option = *options.begin();
  }
  return option;
#else
  return {};
#endif
}

Configuration ActiveHarmony::fetchConfiguration() {
#ifdef AUTOPAS_ENABLE_HARMONY
  const auto traversalOption = fetchTuningParameter(traversalOptionName, _allowedTraversalOptions);
  const auto dataLayoutOption = fetchTuningParameter(dataLayoutOptionName, _allowedDataLayoutOptions);
  const auto newton3Option = fetchTuningParameter(newton3OptionName, _allowedNewton3Options);
  const auto containerOption = *compatibleTraversals::allCompatibleContainers(traversalOption).begin();
  const auto applicableLoadEstimators =
      loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, _allowedLoadEstimatorOptions);
  const auto loadEstimatorOption = fetchTuningParameter(loadEstimatorOptionName, applicableLoadEstimators);

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

  return {containerOption,  cellSizeFactor, traversalOption, loadEstimatorOption,
          dataLayoutOption, newton3Option,  _interactionType};
#else
  return {};
#endif
}

void ActiveHarmony::rejectConfiguration(const Configuration &configuration, bool indefinitely) {
#ifdef AUTOPAS_ENABLE_HARMONY
  // dummy evidence where iteration information is not needed, because we only store the measurement in harmony.
  const Evidence badDummyEvidence{
      0,
      _tuningPhase,
      std::numeric_limits<long>::max(),
  };
  addEvidence(configuration, badDummyEvidence);
  resetHarmony();
#endif
}

bool ActiveHarmony::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                        const EvidenceCollection &evidenceCollection) {
#ifdef AUTOPAS_ENABLE_HARMONY
  // get configurations from server until new configuration with valid newton3 option is found
  bool skipConfig;
  do {
    skipConfig = false;
    if (ah_fetch(htask) < 0) {
      utils::ExceptionHandler::exception("ActiveHarmony::optimizeSuggestions: Error fetching values from server");
    }
    const auto potentialConf = fetchConfiguration();
    // If we already know the performance for this config in this tuning phase skip it.
    if (const auto *potentialConfEvidence = evidenceCollection.getEvidence(potentialConf);
        potentialConfEvidence and not potentialConfEvidence->empty() and
        potentialConfEvidence->back().tuningPhase == _tuningPhase) {
      addEvidence(potentialConf, potentialConfEvidence->back());
      skipConfig = true;
    }

    auto converged = ah_converged(htask);
    if (converged) {
      // set configuration to optimum
      AutoPasLog(TRACE, "ActiveHarmony converged.");
      if (ah_best(htask) != 0) {
        utils::ExceptionHandler::exception("ActiveHarmony::optimizeSuggestions: Error fetching best point.");
      }
      // only accept the configuration if it was in the queue.
      if (std::find(configQueue.begin(), configQueue.end(), potentialConf) == configQueue.end()) {
        skipConfig = true;
      } else {
        return false;
      }
    }

    if (_nonLocalServer) {
      // When using a non-local server, it is possible that only tested configurations are fetched before the search
      // converges.
      // Because this is difficult to test for, the loop is simply ignored for non-local servers.
      return false;
    }
  } while (skipConfig);
#endif
  // ActiveHarmony does no intentional config wipes to stop the tuning phase
  return false;
}

template <class OptionClass>
void ActiveHarmony::configureTuningParameter(hdef_t *hdef, const char *name, const std::set<OptionClass> &options) {
#ifdef AUTOPAS_ENABLE_HARMONY
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
#endif
}

bool ActiveHarmony::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                          const EvidenceCollection &evidenceCollection) {
  resetHarmony();
  // ActiveHarmony does no intentional config wipes to stop the tuning phase
  return false;
}

void ActiveHarmony::resetHarmony() {
#ifdef AUTOPAS_ENABLE_HARMONY
  int rank, commSize;
  AutoPas_MPI_Comm_size(_comm, &commSize);
  AutoPas_MPI_Comm_rank(_comm, &rank);

  if (_mpiDivideAndConquer) {
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
    if (_mpiDivideAndConquer) {
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

      if (_mpiDivideAndConquer) {
        AutoPas_MPI_Barrier(_comm);
      }
    } else {
      // only join a session if using the divideAndConquer mpi strategy, we are not rank 0 and a server is specified.
      if (_mpiDivideAndConquer) {
        AutoPas_MPI_Barrier(_comm);
      }

      // Everybody else may now join the master's new Harmony search.
      htask = ah_join(hdesc, "AutoPas");
      if (htask == nullptr) {
        utils::ExceptionHandler::exception("ActiveHarmony::reset: Error joining session");
      }
    }
  }
#endif
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

void ActiveHarmony::setupTuningParameters(int commSize, hdef_t *hdef) {
#ifdef AUTOPAS_ENABLE_HARMONY
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

  if (_mpiDivideAndConquer) {
    // set the size of the tuning session
    char numbuf[12];
    snprintf(numbuf, sizeof(numbuf), "%d", commSize);
    ah_def_cfg(hdef, "CLIENT_COUNT", numbuf);
  }
#endif
}

bool ActiveHarmony::needsDomainSimilarityStatistics() const { return false; }

TuningStrategyOption ActiveHarmony::getOptionType() const { return TuningStrategyOption::activeHarmony; }
}  // namespace autopas