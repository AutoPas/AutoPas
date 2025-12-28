/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 01.03.2021
 */
#include "Simulation.h"

#include <algorithm>

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/utils/WrapOpenMP.h"

// Declare the main AutoPas class and the computeInteractions() methods with all used functors as extern template
// instantiation. They are instantiated in the respective cpp file inside the templateInstantiations folder.
//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeAutovec *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeAVX *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeSVE *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(ATFunctor *);
#endif
//! @endcond

#include <sys/ioctl.h>
#include <unistd.h>

#include <iostream>

#include "ParticleCommunicator.h"
#include "Thermostat.h"
#include "autopas/utils/MemoryProfiler.h"
#include "autopas/utils/WrapMPI.h"
#include "configuration/MDFlexConfig.h"

namespace {
/**
 * Tries to identify the width of the terminal where the simulation is running.
 * If no width can be identified, the function defaults to 80.
 * @return width of the terminal.
 */
size_t getTerminalWidth() {
  size_t terminalWidth = 0;
  // test all std pipes to get the current terminal width
  for (auto fd : {STDOUT_FILENO, STDIN_FILENO, STDERR_FILENO}) {
    if (isatty(fd)) {
      struct winsize w {};
      ioctl(fd, TIOCGWINSZ, &w);
      terminalWidth = w.ws_col;
      break;
    }
  }

  // if width is still zero try the environment variable COLUMNS
  if (terminalWidth == 0) {
    if (auto *terminalWidthCharArr = std::getenv("COLUMNS")) {
      // this pointer could be used to detect parsing errors via terminalWidthCharArr == end
      // but since we have a fallback further down we are ok if this fails silently.
      char *end{};
      terminalWidth = std::strtol(terminalWidthCharArr, &end, 10);
    }
  }

  // if all of the above fail fall back to a fixed width
  if (terminalWidth == 0) {
    // this seems to be the default width in most terminal windows
    terminalWidth = 80;
  }

  return terminalWidth;
}
}  // namespace

Simulation::Simulation(const MDFlexConfig &configuration,
                       std::shared_ptr<RegularGridDecomposition> &domainDecomposition)
    : _configuration(configuration),
      _domainDecomposition(domainDecomposition),
      _createVtkFiles(not configuration.vtkFileName.value.empty()),
      _vtkWriter(nullptr),
      _totalEnergySensor(configuration.energySensorOption.value) {
  _timers.total.start();
  _timers.initialization.start();
  _totalEnergySensor.startMeasurement();

  // only create the writer if necessary since this also creates the output dir
  if (_createVtkFiles) {
    _vtkWriter =
        std::make_shared<ParallelVtkWriter>(_configuration.vtkFileName.value, _configuration.vtkOutputFolder.value,
                                            std::to_string(_configuration.iterations.value).size());
  }

  const auto rank = _domainDecomposition->getDomainIndex();
  const auto *fillerBeforeSuffix =
      _configuration.outputSuffix.value.empty() or _configuration.outputSuffix.value.front() == '_' ? "" : "_";
  const auto *fillerAfterSuffix =
      _configuration.outputSuffix.value.empty() or _configuration.outputSuffix.value.back() == '_' ? "" : "_";
  const auto outputSuffix =
      "Rank" + std::to_string(rank) + fillerBeforeSuffix + _configuration.outputSuffix.value + fillerAfterSuffix;

  if (rank == 0) {
    _globalLogger = std::make_unique<GlobalVariableLogger>(outputSuffix);
  }

  if (_configuration.logFileName.value.empty()) {
    _outputStream = &std::cout;
  } else {
    _logFile = std::make_shared<std::ofstream>();
    _logFile->open(_configuration.logFileName.value + "_" + outputSuffix);
    _outputStream = &(*_logFile);
  }

  _autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(*_outputStream);
  _autoPasContainer->setAllowedCellSizeFactors(*_configuration.cellSizeFactors.value);
  _autoPasContainer->setAllowedContainers(_configuration.containerOptions.value);

  if (_configuration.getInteractionTypes().empty()) {
    std::string functorName{};
    std::tie(_configuration.functorOption.value, functorName) =
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
        std::make_pair(MDFlexConfig::FunctorOption::lj12_6_AVX, "Lennard-Jones AVX Functor.");
#elif defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
        std::make_pair(MDFlexConfig::FunctorOption::lj12_6_SVE, "Lennard-Jones SVE Functor.");
#else
        std::make_pair(MDFlexConfig::FunctorOption::lj12_6, "Lennard-Jones AutoVec Functor.");
#endif
    _configuration.addInteractionType(autopas::InteractionTypeOption::pairwise);
    std::cout << "WARNING: No functor was specified. Defaulting to " << functorName << std::endl;
  }

  _autoPasContainer->setAllowedInteractionTypeOptions(_configuration.getInteractionTypes());

  // Pairwise specific options
  _autoPasContainer->setAllowedDataLayouts(_configuration.dataLayoutOptions.value,
                                           autopas::InteractionTypeOption::pairwise);
  _autoPasContainer->setAllowedNewton3Options(_configuration.newton3Options.value,
                                              autopas::InteractionTypeOption::pairwise);
  _autoPasContainer->setAllowedTraversals(_configuration.traversalOptions.value,
                                          autopas::InteractionTypeOption::pairwise);
  _autoPasContainer->setAllowedLoadEstimators(_configuration.loadEstimatorOptions.value);
  // Triwise specific options
  _autoPasContainer->setAllowedDataLayouts(_configuration.dataLayoutOptions3B.value,
                                           autopas::InteractionTypeOption::triwise);
  _autoPasContainer->setAllowedNewton3Options(_configuration.newton3Options3B.value,
                                              autopas::InteractionTypeOption::triwise);
  _autoPasContainer->setAllowedTraversals(_configuration.traversalOptions3B.value,
                                          autopas::InteractionTypeOption::triwise);
  // General options
  _autoPasContainer->setBoxMin(_domainDecomposition->getLocalBoxMin());
  _autoPasContainer->setBoxMax(_domainDecomposition->getLocalBoxMax());
  _autoPasContainer->setCutoff(_configuration.cutoff.value);
  _autoPasContainer->setRelativeOptimumRange(_configuration.relativeOptimumRange.value);
  _autoPasContainer->setMaxTuningPhasesWithoutTest(_configuration.maxTuningPhasesWithoutTest.value);
  _autoPasContainer->setRelativeBlacklistRange(_configuration.relativeBlacklistRange.value);
  _autoPasContainer->setEvidenceFirstPrediction(_configuration.evidenceFirstPrediction.value);
  _autoPasContainer->setExtrapolationMethodOption(_configuration.extrapolationMethodOption.value);
  _autoPasContainer->setNumSamples(_configuration.tuningSamples.value);
  _autoPasContainer->setEarlyStoppingFactor(_configuration.earlyStoppingFactor.value);
  _autoPasContainer->setMaxEvidence(_configuration.tuningMaxEvidence.value);
  _autoPasContainer->setRuleFileName(_configuration.ruleFilename.value);
  _autoPasContainer->setFuzzyRuleFileName(_configuration.fuzzyRuleFilename.value);
  _autoPasContainer->setSelectorStrategy(_configuration.selectorStrategy.value);
  _autoPasContainer->setTuningInterval(_configuration.tuningInterval.value);
  _autoPasContainer->setTuningStrategyOption(_configuration.tuningStrategyOptions.value);
  _autoPasContainer->setTuningMetricOption(_configuration.tuningMetricOption.value);
  _autoPasContainer->setUseLOESSSmoothening(_configuration.useLOESSSmoothening.value);
  _autoPasContainer->setEnergySensorOption(_configuration.energySensorOption.value);
  _autoPasContainer->setMPITuningMaxDifferenceForBucket(_configuration.MPITuningMaxDifferenceForBucket.value);
  _autoPasContainer->setMPITuningWeightForMaxDensity(_configuration.MPITuningWeightForMaxDensity.value);
  _autoPasContainer->setVerletClusterSize(_configuration.verletClusterSize.value);
  _autoPasContainer->setVerletRebuildFrequency(_configuration.verletRebuildFrequency.value);
  _autoPasContainer->setVerletSkin(_configuration.verletSkinRadius.value);
  _autoPasContainer->setAcquisitionFunction(_configuration.acquisitionFunctionOption.value);
  _autoPasContainer->setUseTuningLogger(_configuration.useTuningLogger.value);
  _autoPasContainer->setSortingThreshold(_configuration.sortingThreshold.value);
  _autoPasContainer->setOutputSuffix(outputSuffix);
  autopas::Logger::get()->set_level(_configuration.logLevel.value);

  _autoPasContainer->init();

  // Throw an error if there is not more than one configuration to test in the search space but more than one tuning
  // phase is requested
  if (_autoPasContainer->searchSpaceIsTrivial() and _configuration.tuningPhases.value > 0) {
    throw std::runtime_error(
        "Search space must not be trivial if the simulation time is limited by the number tuning phases");
  }

  // Load particles from the config file (object generators, checkpoint)
  loadParticles();

  if (_configuration.useThermostat.value and _configuration.deltaT.value != 0) {
    if (_configuration.addBrownianMotion.value) {
      Thermostat::addBrownianMotion(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                    _configuration.initTemperature.value);
    }

    // Set the simulation directly to the desired initial temperature.
    Thermostat::apply(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                      _configuration.initTemperature.value, std::numeric_limits<double>::max());
  }

  _timers.initialization.stop();
}

void Simulation::finalize() {
  _timers.total.stop();
  _totalEnergySensor.endMeasurement();
  autopas::AutoPas_MPI_Barrier(AUTOPAS_MPI_COMM_WORLD);

  logSimulationState();
  logMeasurements();
}

void Simulation::run() {
  _timers.simulate.start();
  while (needsMoreIterations()) {
    if (_createVtkFiles and _iteration % _configuration.vtkWriteFrequency.value == 0) {
      _timers.vtk.start();
      _vtkWriter->recordTimestep(_iteration, *_autoPasContainer, *_domainDecomposition);
      _timers.vtk.stop();
    }

    _timers.computationalLoad.start();
    if (_configuration.deltaT.value != 0 and not _simulationIsPaused) {
      updatePositionsAndResetForces();
#if MD_FLEXIBLE_MODE == MULTISITE
      updateQuaternions();
#endif
    }

    // We update the container, even if dt=0, to bump the iteration counter, which is needed to ensure containers can
    // still be rebuilt in frozen scenarios e.g. for algorithm performance data gathering purposes. Also, it bumps the
    // iteration counter which can be used to uniquely identify functor calls.
    _timers.updateContainer.start();
    auto emigrants = _autoPasContainer->updateContainer();
    _timers.updateContainer.stop();

    if (_configuration.deltaT.value != 0 and not _simulationIsPaused) {
      const auto computationalLoad = static_cast<double>(_timers.computationalLoad.stop());

      // periodically resize box for MPI load balancing
      if (_iteration % _configuration.loadBalancingInterval.value == 0) {
        _timers.loadBalancing.start();
        _domainDecomposition->update(computationalLoad);
        auto additionalEmigrants = _autoPasContainer->resizeBox(_domainDecomposition->getLocalBoxMin(),
                                                                _domainDecomposition->getLocalBoxMax());
        // If the boundaries shifted, particles that were thrown out by updateContainer() previously might now be in the
        // container again.
        // Reinsert emigrants if they are now inside the domain and mark local copies as dummy,
        // so that remove_if can erase them after.
        const auto &boxMin = _autoPasContainer->getBoxMin();
        const auto &boxMax = _autoPasContainer->getBoxMax();
        _autoPasContainer->addParticlesIf(emigrants, [&](auto &p) {
          if (autopas::utils::inBox(p.getR(), boxMin, boxMax)) {
            // This only changes the ownership state in the emigrants vector, not in AutoPas
            p.setOwnershipState(autopas::OwnershipState::dummy);
            return true;
          }
          return false;
        });

        emigrants.erase(std::remove_if(emigrants.begin(), emigrants.end(), [&](const auto &p) { return p.isDummy(); }),
                        emigrants.end());

        emigrants.insert(emigrants.end(), additionalEmigrants.begin(), additionalEmigrants.end());
        _timers.loadBalancing.stop();
      }

      _timers.migratingParticleExchange.start();
      _domainDecomposition->exchangeMigratingParticles(*_autoPasContainer, emigrants);
      _timers.migratingParticleExchange.stop();

      _timers.reflectParticlesAtBoundaries.start();
      _domainDecomposition->reflectParticlesAtBoundaries(*_autoPasContainer,
                                                         *_configuration.getParticlePropertiesLibrary());
      _timers.reflectParticlesAtBoundaries.stop();

      _timers.haloParticleExchange.start();
      _domainDecomposition->exchangeHaloParticles(*_autoPasContainer);
      _timers.haloParticleExchange.stop();

      _timers.computationalLoad.start();
    }

    updateInteractionForces();

    if (_configuration.pauseSimulationDuringTuning.value) {
      // If PauseSimulationDuringTuning is enabled we need to update the _simulationIsPaused flag
      updateSimulationPauseState();
    }

    if (_configuration.deltaT.value != 0 and not _simulationIsPaused) {
      updateVelocities();
#if MD_FLEXIBLE_MODE == MULTISITE
      updateAngularVelocities();
#endif
      updateThermostat();
    }
    _timers.computationalLoad.stop();
#ifdef MD_FLEXIBLE_CALC_GLOBALS
    // Summing the potential energy over all MPI ranks
    double potentialEnergyOverMPIRanks{}, virialSumOverMPIRanks{};
    autopas::AutoPas_MPI_Reduce(&_totalPotentialEnergy, &potentialEnergyOverMPIRanks, 1, AUTOPAS_MPI_DOUBLE,
                                AUTOPAS_MPI_SUM, 0, AUTOPAS_MPI_COMM_WORLD);
    autopas::AutoPas_MPI_Reduce(&_totalVirialSum, &virialSumOverMPIRanks, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM, 0,
                                AUTOPAS_MPI_COMM_WORLD);
    if (_domainDecomposition->getDomainIndex() == 0) {
      _globalLogger->logGlobals(_iteration, potentialEnergyOverMPIRanks, virialSumOverMPIRanks);
    }
    _totalPotentialEnergy = 0.;
    _totalVirialSum = 0.;
#endif
    if (not _simulationIsPaused) {
      ++_iteration;
    }

    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage on rank " << _domainDecomposition->getDomainIndex() << ": "
                << autopas::memoryProfiler::currentMemoryUsage() << " kB\n";
    }

    if (_domainDecomposition->getDomainIndex() == 0) {
      auto [maxIterationsEstimate, maxIterationsIsPrecise] = estimateNumberOfIterations();
      if (not _configuration.dontShowProgressBar.value) {
        printProgress(_iteration, maxIterationsEstimate, maxIterationsIsPrecise);
      }
    }
  }
  _timers.simulate.stop();

  // Record last state of simulation.
  if (_createVtkFiles) {
    _vtkWriter->recordTimestep(_iteration, *_autoPasContainer, *_domainDecomposition);
  }
}

std::tuple<size_t, bool> Simulation::estimateNumberOfIterations() const {
  if (_configuration.tuningPhases.value > 0) {
    const size_t configsTestedPerTuningPhase = [&]() {
      if (std::any_of(_configuration.tuningStrategyOptions.value.begin(),
                      _configuration.tuningStrategyOptions.value.end(), [](const auto &stratOpt) {
                        return stratOpt == autopas::TuningStrategyOption::bayesianSearch or
                               stratOpt == autopas::TuningStrategyOption::bayesianClusterSearch or
                               stratOpt == autopas::TuningStrategyOption::randomSearch;
                      })) {
        return static_cast<size_t>(_configuration.tuningMaxEvidence.value);
      } else {
        // @TODO: this can be improved by considering the tuning strategy
        //      or averaging number of iterations per tuning phase and dynamically adapt prediction

        // This estimate is only valid for full search and no restrictions on the cartesian product.
        // add static to only evaluate this once
        const size_t searchSpaceSizePairwise =
            _configuration.getInteractionTypes().count(autopas::InteractionTypeOption::pairwise) == 0
                ? 0
                : autopas::SearchSpaceGenerators::cartesianProduct(
                      _configuration.containerOptions.value, _configuration.traversalOptions.value,
                      _configuration.loadEstimatorOptions.value, _configuration.dataLayoutOptions.value,
                      _configuration.newton3Options.value, _configuration.cellSizeFactors.value.get(),
                      autopas::InteractionTypeOption::pairwise)
                      .size();

        const size_t searchSpaceSizeTriwise =
            _configuration.getInteractionTypes().count(autopas::InteractionTypeOption::triwise) == 0
                ? 0
                : autopas::SearchSpaceGenerators::cartesianProduct(
                      _configuration.containerOptions.value, _configuration.traversalOptions3B.value,
                      _configuration.loadEstimatorOptions.value, _configuration.dataLayoutOptions3B.value,
                      _configuration.newton3Options3B.value, _configuration.cellSizeFactors.value.get(),
                      autopas::InteractionTypeOption::triwise)
                      .size();

        return std::max(searchSpaceSizePairwise, searchSpaceSizeTriwise);
      }
    }();
    // non-tuning iterations + tuning iterations + one iteration after last phase
    auto estimate =
        (_configuration.tuningPhases.value - 1) * _configuration.tuningInterval.value +
        (_configuration.tuningPhases.value * _configuration.tuningSamples.value * configsTestedPerTuningPhase) + 1;
    return {estimate, false};
  } else {
    return {_configuration.iterations.value, true};
  }
}

void Simulation::printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise) {
  // percentage of iterations complete
  const double fractionDone = static_cast<double>(iterationProgress) / static_cast<double>(maxIterations);

  // length of the number of maxIterations
  const auto numCharsOfMaxIterations = static_cast<int>(std::to_string(maxIterations).size());

  // trailing information string
  std::stringstream info;
  info << std::setw(3) << std::round(fractionDone * 100) << "% " << std::setw(numCharsOfMaxIterations)
       << iterationProgress << "/";
  if (not maxIsPrecise) {
    info << "~";
  }
  info << maxIterations;

  // actual progress bar
  std::stringstream progressbar;
  progressbar << "[";
  // get current terminal width
  const auto terminalWidth = getTerminalWidth();

  // the bar should fill the terminal window so subtract everything else (-2 for "] ")
  const int maxBarWidth = static_cast<int>(terminalWidth - info.str().size() - progressbar.str().size() - 2ul);
  // sanity check for underflow
  if (maxBarWidth > terminalWidth) {
    std::cerr << "Warning! Terminal width appears to be too small or could not be read. Disabling progress bar."
              << std::endl;
    _configuration.dontShowProgressBar.value = true;
    return;
  }
  const auto barWidth =
      std::max(std::min(static_cast<decltype(maxBarWidth)>(maxBarWidth * (fractionDone)), maxBarWidth), 1);
  // don't print arrow tip if >= 100%
  if (iterationProgress >= maxIterations) {
    progressbar << std::string(barWidth, '=');
  } else {
    progressbar << std::string(barWidth - 1, '=') << '>' << std::string(maxBarWidth - barWidth, ' ');
  }
  progressbar << "] ";
  // clear current line (=delete previous progress bar)
  std::cout << std::string(terminalWidth, '\r');
  // print everything
  std::cout << progressbar.str() << info.str() << std::flush;
}

std::string Simulation::timerToString(const std::string &name, long timeNS, int numberWidth, long maxTime) {
  // only print timers that were actually used
  if (timeNS == 0) {
    return "";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(_floatStringPrecision) << name << " : " << std::setw(numberWidth) << std::right
     << timeNS
     << " ns ("
     // min width of the representation of seconds is numberWidth - 9 (from conversion) + 4 (for dot and digits after)
     << std::setw(numberWidth - 5) << ((double)timeNS * 1e-9) << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeNS / (double)maxTime * 100) << "%";
  }
  ss << "\n";
  return ss.str();
}

void Simulation::updatePositionsAndResetForces() {
  _timers.positionUpdate.start();
  TimeDiscretization::calculatePositionsAndResetForces(
      *_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()), _configuration.deltaT.value,
      _configuration.globalForce.value, _configuration.fastParticlesThrow.value);
  _timers.positionUpdate.stop();
}

void Simulation::updateQuaternions() {
  _timers.quaternionUpdate.start();
  TimeDiscretization::calculateQuaternionsAndResetTorques(
      *_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()), _configuration.deltaT.value,
      _configuration.globalForce.value);
  _timers.quaternionUpdate.stop();
}

void Simulation::updateInteractionForces() {
  _timers.forceUpdateTotal.start();

  _previousIterationWasTuningIteration = _currentIterationIsTuningIteration;
  _currentIterationIsTuningIteration = false;
  long timeIteration = 0;

  // Calculate pairwise forces
  if (_configuration.getInteractionTypes().count(autopas::InteractionTypeOption::pairwise)) {
    _timers.forceUpdatePairwise.start();
    _currentIterationIsTuningIteration |= calculatePairwiseForces();
    timeIteration += _timers.forceUpdatePairwise.stop();
  }
  // Calculate triwise forces
  if (_configuration.getInteractionTypes().count(autopas::InteractionTypeOption::triwise)) {
    _timers.forceUpdateTriwise.start();
    _currentIterationIsTuningIteration |= calculateTriwiseForces();
    timeIteration += _timers.forceUpdateTriwise.stop();
  }

  // count time spent for tuning
  if (_currentIterationIsTuningIteration) {
    _timers.forceUpdateTuning.addTime(timeIteration);
    ++_numTuningIterations;
  } else {
    _timers.forceUpdateNonTuning.addTime(timeIteration);
    // if the previous iteration was a tuning iteration and the current one is not
    // we have reached the end of a tuning phase
    if (_previousIterationWasTuningIteration) {
      ++_numTuningPhasesCompleted;
    }
  }

  _timers.forceUpdateTotal.stop();
}

void Simulation::updateVelocities() {
  const double deltaT = _configuration.deltaT.value;

  if (deltaT != 0) {
    _timers.velocityUpdate.start();
    TimeDiscretization::calculateVelocities(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                            deltaT);
    _timers.velocityUpdate.stop();
  }
}

void Simulation::updateAngularVelocities() {
  const double deltaT = _configuration.deltaT.value;

  _timers.angularVelocityUpdate.start();
  TimeDiscretization::calculateAngularVelocities(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                                 deltaT);
  _timers.angularVelocityUpdate.stop();
}

void Simulation::updateThermostat() {
  if (_configuration.useThermostat.value and (_iteration % _configuration.thermostatInterval.value) == 0) {
    _timers.thermostat.start();
    Thermostat::apply(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                      _configuration.targetTemperature.value, _configuration.deltaTemp.value);
    _timers.thermostat.stop();
  }
}

long Simulation::accumulateTime(const long &time) {
  long reducedTime{};
  autopas::AutoPas_MPI_Reduce(&time, &reducedTime, 1, AUTOPAS_MPI_LONG, AUTOPAS_MPI_SUM, 0, AUTOPAS_MPI_COMM_WORLD);

  return reducedTime;
}

bool Simulation::calculatePairwiseForces() {
  const auto wasTuningIteration = applyWithChosenFunctor<bool>([&](auto &&functor) {
    auto isTuningIteration = _autoPasContainer->computeInteractions(&functor);
#ifdef MD_FLEXIBLE_CALC_GLOBALS
    _totalPotentialEnergy += functor.getPotentialEnergy();
    _totalVirialSum += functor.getVirial();
#endif
    return isTuningIteration;
  });
  return wasTuningIteration;
}

bool Simulation::calculateTriwiseForces() {
  const auto wasTuningIteration = applyWithChosenFunctor3B<bool>([&](auto &&functor) {
    auto isTuningIteration = _autoPasContainer->computeInteractions(&functor);
#ifdef MD_FLEXIBLE_CALC_GLOBALS
    _totalPotentialEnergy += functor.getPotentialEnergy();
    _totalVirialSum += functor.getVirial();
#endif
    return isTuningIteration;
  });
  return wasTuningIteration;
}

void Simulation::calculateGlobalForces(const std::array<double, 3> &globalForce) {
  AUTOPAS_OPENMP(parallel shared(_autoPasContainer))
  for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    particle->addF(globalForce);
  }
}

void Simulation::logSimulationState() {
  size_t totalNumberOfParticles{0ul}, ownedParticles{0ul}, haloParticles{0ul};

  size_t particleCount = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo);
  autopas::AutoPas_MPI_Allreduce(&particleCount, &totalNumberOfParticles, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);

  particleCount = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned);
  autopas::AutoPas_MPI_Allreduce(&particleCount, &ownedParticles, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);

  particleCount = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo);
  autopas::AutoPas_MPI_Allreduce(&particleCount, &haloParticles, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);

  if (_domainDecomposition->getDomainIndex() == 0) {
    std::cout << "\n\n"
              << "Total number of particles at the end of Simulation: " << totalNumberOfParticles << "\n"
              << "Owned: " << ownedParticles << "\n"
              << "Halo : " << haloParticles << "\n";
  }
}

void Simulation::updateSimulationPauseState() {
  // If we are at the beginning of a tuning phase, we need to freeze the simulation
  if (_currentIterationIsTuningIteration and (not _previousIterationWasTuningIteration)) {
    std::cout << "Iteration " << _iteration << ": Freezing simulation for tuning phase. Starting tuning phase...\n";
    _simulationIsPaused = true;
  }

  // If we are at the end of a tuning phase, we need to resume the simulation
  if (_previousIterationWasTuningIteration and (not _currentIterationIsTuningIteration)) {
    std::cout << "Iteration " << _iteration << ": Resuming simulation after tuning phase.\n";

    // reset the forces which accumulated during the tuning phase
    for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      particle->setF(_configuration.globalForce.value);
    }

    // calculate the forces of the latest iteration again
    updateInteractionForces();

    _simulationIsPaused = false;
  }
}

void Simulation::logMeasurements() {
  const long positionUpdate = accumulateTime(_timers.positionUpdate.getTotalTime());
  const long quaternionUpdate = accumulateTime(_timers.quaternionUpdate.getTotalTime());
  const long updateContainer = accumulateTime(_timers.updateContainer.getTotalTime());
  const long forceUpdateTotal = accumulateTime(_timers.forceUpdateTotal.getTotalTime());
  const long forceUpdatePairwise = accumulateTime(_timers.forceUpdatePairwise.getTotalTime());
  const long forceUpdateTriwise = accumulateTime(_timers.forceUpdateTriwise.getTotalTime());
  const long forceUpdateTuning = accumulateTime(_timers.forceUpdateTuning.getTotalTime());
  const long forceUpdateNonTuning = accumulateTime(_timers.forceUpdateNonTuning.getTotalTime());
  const long velocityUpdate = accumulateTime(_timers.velocityUpdate.getTotalTime());
  const long angularVelocityUpdate = accumulateTime(_timers.angularVelocityUpdate.getTotalTime());
  const long simulate = accumulateTime(_timers.simulate.getTotalTime());
  const long vtk = accumulateTime(_timers.vtk.getTotalTime());
  const long initialization = accumulateTime(_timers.initialization.getTotalTime());
  const long total = accumulateTime(_timers.total.getTotalTime());
  const long thermostat = accumulateTime(_timers.thermostat.getTotalTime());
  const long haloParticleExchange = accumulateTime(_timers.haloParticleExchange.getTotalTime());
  const long reflectParticlesAtBoundaries = accumulateTime(_timers.reflectParticlesAtBoundaries.getTotalTime());
  const long migratingParticleExchange = accumulateTime(_timers.migratingParticleExchange.getTotalTime());
  const long loadBalancing = accumulateTime(_timers.loadBalancing.getTotalTime());
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
  double totalEnergy = _totalEnergySensor.getJoules();
  autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &totalEnergy, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);
#endif

  if (_domainDecomposition->getDomainIndex() == 0) {
    const long wallClockTime = _timers.total.getTotalTime();
    // use the two potentially largest timers to determine the number of chars needed
    const auto maximumNumberOfDigits =
        static_cast<int>(std::max(std::to_string(total).length(), std::to_string(wallClockTime).length()));
    std::cout << "Measurements:\n";
    std::cout << timerToString("Total accumulated                 ", total, maximumNumberOfDigits);
    std::cout << timerToString("  Initialization                  ", initialization, maximumNumberOfDigits, total);
    std::cout << timerToString("  Simulate                        ", simulate, maximumNumberOfDigits, total);
    std::cout << timerToString("    PositionUpdate                ", positionUpdate, maximumNumberOfDigits, simulate);
#if MD_FLEXIBLE_MODE == MULTISITE
    std::cout << timerToString("    QuaternionUpdate              ", quaternionUpdate, maximumNumberOfDigits, simulate);
#endif
    std::cout << timerToString("    UpdateContainer               ", updateContainer, maximumNumberOfDigits, simulate);
    std::cout << timerToString("    Boundaries                    ", haloParticleExchange + migratingParticleExchange,
                               maximumNumberOfDigits, simulate);
    std::cout << timerToString("      HaloParticleExchange        ", haloParticleExchange, maximumNumberOfDigits,
                               haloParticleExchange + reflectParticlesAtBoundaries + migratingParticleExchange);
    std::cout << timerToString("      ReflectParticlesAtBoundaries", reflectParticlesAtBoundaries,
                               maximumNumberOfDigits,
                               haloParticleExchange + reflectParticlesAtBoundaries + migratingParticleExchange);
    std::cout << timerToString("      MigratingParticleExchange   ", migratingParticleExchange, maximumNumberOfDigits,
                               haloParticleExchange + reflectParticlesAtBoundaries + migratingParticleExchange);
    std::cout << timerToString("    ForceUpdateTotal              ", forceUpdateTotal, maximumNumberOfDigits, simulate);
    std::cout << timerToString("      Tuning                      ", forceUpdateTuning, maximumNumberOfDigits,
                               forceUpdateTotal);
    std::cout << timerToString("      ForceUpdateTuning           ", forceUpdateTuning, maximumNumberOfDigits,
                               forceUpdateTotal);
    std::cout << timerToString("      ForceUpdateNonTuning        ", forceUpdateNonTuning, maximumNumberOfDigits,
                               forceUpdateTotal);
    std::cout << timerToString("    VelocityUpdate                ", velocityUpdate, maximumNumberOfDigits, simulate);
#if MD_FLEXIBLE_MODE == MULTISITE
    std::cout << timerToString("    AngularVelocityUpdate         ", angularVelocityUpdate, maximumNumberOfDigits,
                               simulate);
#endif
    std::cout << timerToString("    Thermostat                    ", thermostat, maximumNumberOfDigits, simulate);
    std::cout << timerToString("    Vtk                           ", vtk, maximumNumberOfDigits, simulate);
    std::cout << timerToString("    LoadBalancing                 ", loadBalancing, maximumNumberOfDigits, simulate);
    std::cout << timerToString("  One iteration                 ", simulate / static_cast<long>(_iteration),
                               maximumNumberOfDigits, total);

    std::cout << timerToString("Total wall-clock time             ", wallClockTime, maximumNumberOfDigits, total);
    std::cout << "\n";

    std::cout << "Tuning iterations                  : " << _numTuningIterations << " / " << _iteration << " = "
              << (static_cast<double>(_numTuningIterations) / static_cast<double>(_iteration) * 100.) << "%"
              << "\n";

    auto mfups =
        static_cast<double>(_autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned) * _iteration) *
        1e-6 / (static_cast<double>(forceUpdateTotal) * 1e-9);  // 1e-9 for ns to s, 1e-6 for M in MFUPs
    std::cout << "MFUPs/sec                          : " << mfups << "\n";
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    std::cout << "Mean Rebuild Frequency               : " << _autoPasContainer->getMeanRebuildFrequency() << "\n";
#endif
#ifdef AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
    std::cout << "Total Energy Consumed (in Joules)    : " << totalEnergy << "\n";
#endif
  }
}

bool Simulation::needsMoreIterations() const {
  return _iteration < _configuration.iterations.value or _numTuningPhasesCompleted < _configuration.tuningPhases.value;
}

void Simulation::checkNumParticles(size_t expectedNumParticlesGlobal, size_t numParticlesCurrentlyMigratingLocal,
                                   int lineNumber) {
  if (std::all_of(_configuration.boundaryOption.value.begin(), _configuration.boundaryOption.value.end(),
                  [](const auto &boundary) {
                    return boundary == options::BoundaryTypeOption::periodic or
                           boundary == options::BoundaryTypeOption::reflective;
                  })) {
    const auto numParticlesNowLocal = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned);
    std::array<size_t, 2> sendBuffer{numParticlesNowLocal, numParticlesCurrentlyMigratingLocal};
    std::array<size_t, sendBuffer.size()> receiveBuffer{};
    autopas::AutoPas_MPI_Reduce(sendBuffer.data(), receiveBuffer.data(), receiveBuffer.size(),
                                AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0, AUTOPAS_MPI_COMM_WORLD);
    const auto &[numParticlesNowTotal, numParticlesMigratingTotal] = receiveBuffer;
    if (expectedNumParticlesGlobal != numParticlesNowTotal + numParticlesMigratingTotal) {
      const auto myRank = _domainDecomposition->getDomainIndex();
      std::stringstream ss;
      // clang-format off
      ss << "Rank " << myRank << " Line " << lineNumber
         << ": Particles Lost! All Boundaries are periodic but the number of particles changed:"
         << "Expected        : " << expectedNumParticlesGlobal
         << "Actual          : " << (numParticlesNowTotal + numParticlesMigratingTotal)
         << "  in containers : " << numParticlesNowTotal
         << "  migrating     : " << numParticlesMigratingTotal
         << std::endl;
      // clang-format on
      throw std::runtime_error(ss.str());
    }
  }
}

void Simulation::loadParticles() {
  // Store how many particles are in this config object before removing them.
  const auto numParticlesInConfigLocally = _configuration.particles.size();
  // When loading checkpoints, the config file might contain particles that do not belong to this rank,
  // because rank boundaries don't align with those at the end of the previous simulation due to dynamic load balancing.
  // Idea: Only load what belongs here and send the rest away.
  _autoPasContainer->addParticlesIf(_configuration.particles, [&](auto &p) {
    if (_domainDecomposition->isInsideLocalDomain(p.getR())) {
      // Mark particle in vector as dummy, so we know it has been inserted.
      p.setOwnershipState(autopas::OwnershipState::dummy);
      return true;
    }
    return false;
  });

  // Remove what has been inserted. Everything that remains does not belong into this rank.
  _configuration.particles.erase(std::remove_if(_configuration.particles.begin(), _configuration.particles.end(),
                                                [&](const auto &p) { return p.isDummy(); }),
                                 _configuration.particles.end());

  // Send all remaining particles to all ranks
  // TODO: This is not optimal but since this only happens once upon initialization it is not too bad.
  //       Nevertheless it could be improved by determining which particle has to go to which rank.
  const auto rank = _domainDecomposition->getDomainIndex();
  ParticleCommunicator particleCommunicator(_domainDecomposition->getCommunicator());
  for (int receiverRank = 0; receiverRank < _domainDecomposition->getNumberOfSubdomains(); ++receiverRank) {
    // don't send to ourselves
    if (receiverRank == rank) {
      continue;
    }
    particleCommunicator.sendParticles(_configuration.particles, receiverRank);
  }
  // Erase all locally stored particles. They don't belong to this rank and have been sent away.
  _configuration.flushParticles();

  // Receive particles from all other ranks.
  for (int senderRank = 0; senderRank < _domainDecomposition->getNumberOfSubdomains(); ++senderRank) {
    // don't send to ourselves
    if (senderRank == rank) {
      continue;
    }
    particleCommunicator.receiveParticles(_configuration.particles, senderRank);
  }
  particleCommunicator.waitForSendRequests();

  // Add all new particles that belong in this rank
  _autoPasContainer->addParticlesIf(_configuration.particles,
                                    [&](auto &p) { return _domainDecomposition->isInsideLocalDomain(p.getR()); });
  // cleanup
  _configuration.flushParticles();

  // Output and sanity checks
  const size_t numParticlesLocally = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned);
  std::cout << "Number of particles at initialization "
            // align outputs based on the max number of ranks
            << "on rank " << std::setw(std::to_string(_domainDecomposition->getNumberOfSubdomains()).length())
            << std::right << rank << ": " << numParticlesLocally << "\n";

  // Local unnamed struct to pack data for MPI
  struct {
    size_t numParticlesAdded;
    size_t numParticlesInConfig;
  } dataPackage{numParticlesLocally, numParticlesInConfigLocally};
  // Let rank 0 also report the global number of particles
  if (rank == 0) {
    autopas::AutoPas_MPI_Reduce(AUTOPAS_MPI_IN_PLACE, &dataPackage, 2, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0,
                                _domainDecomposition->getCommunicator());
    std::cout << "Number of particles at initialization globally"
              // align ":" with the messages above
              << std::setw(std::to_string(_domainDecomposition->getNumberOfSubdomains()).length()) << ""
              << ": " << dataPackage.numParticlesAdded << "\n";
    // Sanity check that on a global scope all particles have been loaded
    if (dataPackage.numParticlesAdded != dataPackage.numParticlesInConfig) {
      throw std::runtime_error(
          "Simulation::loadParticles(): "
          "Not all particles from the configuration file could be added to AutoPas!\n"
          "Configuration : " +
          std::to_string(dataPackage.numParticlesInConfig) +
          "\n"
          "Added globally: " +
          std::to_string(dataPackage.numParticlesAdded));
    }
  } else {
    // In-place reduce needs different calls on root vs rest...
    autopas::AutoPas_MPI_Reduce(&dataPackage, nullptr, 2, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0,
                                _domainDecomposition->getCommunicator());
  }
}

template <class ReturnType, class FunctionType>
ReturnType Simulation::applyWithChosenFunctor(FunctionType f) {
  const double cutoff = _configuration.cutoff.value;
  auto &particlePropertiesLibrary = *_configuration.getParticlePropertiesLibrary();
  switch (_configuration.functorOption.value) {
    case MDFlexConfig::FunctorOption::lj12_6: {
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
      return f(LJFunctorTypeAutovec{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor AutoVec. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AUTOVEC=ON`.");
#endif
    }
    case MDFlexConfig::FunctorOption::lj12_6_AVX: {
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
      return f(LJFunctorTypeAVX{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor AVX. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AVX=ON`.");
#endif
    }
    case MDFlexConfig::FunctorOption::lj12_6_SVE: {
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
      return f(LJFunctorTypeSVE{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor SVE. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_SVE=ON`.");
#endif
    }
    default: {
      throw std::runtime_error("Unknown pairwise functor choice" +
                               std::to_string(static_cast<int>(_configuration.functorOption.value)));
    }
  }
}

template <class ReturnType, class FunctionType>
ReturnType Simulation::applyWithChosenFunctor3B(FunctionType f) {
  const double cutoff = _configuration.cutoff.value;
  auto &particlePropertiesLibrary = *_configuration.getParticlePropertiesLibrary();
  switch (_configuration.functorOption3B.value) {
    case MDFlexConfig::FunctorOption3B::at: {
#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC)
      return f(ATMFunctor{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for AxilrodTellerMuto Functor. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AT_AUTOVEC=ON`.");
#endif
    }
    default: {
      throw std::runtime_error("Unknown triwise functor choice" +
                               std::to_string(static_cast<int>(_configuration.functorOption3B.value)));
    }
  }
}
