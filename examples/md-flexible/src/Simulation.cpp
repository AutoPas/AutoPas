/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 01.03.2021
 */
#include "Simulation.h"

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"

// Declare the main AutoPas class and the iteratePairwise() methods with all used functors as extern template
// instantiation. They are instantiated in the respective cpp file inside the templateInstantiations folder.
//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
extern template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::LJFunctor<ParticleType, true, true> *);
extern template bool autopas::AutoPas<ParticleType>::iteratePairwise(
    autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true> *);
extern template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::LJFunctorAVX<ParticleType, true, true> *);
extern template bool autopas::AutoPas<ParticleType>::iteratePairwise(autopas::FlopCounterFunctor<ParticleType> *);
//! @endcond

#include <sys/ioctl.h>
#include <unistd.h>

#include <fstream>
#include <iostream>

#include "BoundaryConditions.h"
#include "Thermostat.h"
#include "TimeDiscretization.h"
#include "configuration/MDFlexConfig.h"
#include "src/ParticleSerializationTools.h"

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
    if (auto *teminalWidthCharArr = std::getenv("COLUMNS")) {
      terminalWidth = atoi(teminalWidthCharArr);
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

Simulation::Simulation(const MDFlexConfig &configuration, RegularGridDecomposition &domainDecomposition)
    : _configuration(configuration),
      _domainDecomposition(domainDecomposition),
      _createVtkFiles(not configuration.vtkFileName.value.empty()),
      _vtkWriter(std::make_shared<ParallelVtkWriter>(_configuration.vtkFileName.value,
                                                     _configuration.vtkOutputFolder.value,
                                                     std::to_string(_configuration.iterations.value).size())) {
  _timers.total.start();
  _timers.initialization.start();

  if (_configuration.logFileName.value.empty()) {
    _outputStream = &std::cout;
  } else {
    _logFile = std::make_shared<std::ofstream>();
    _logFile->open(_configuration.logFileName.value);
    _outputStream = &(*_logFile);
  }

  _autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(*_outputStream);
  _autoPasContainer->setAllowedCellSizeFactors(*_configuration.cellSizeFactors.value);
  _autoPasContainer->setAllowedContainers(_configuration.containerOptions.value);
  _autoPasContainer->setAllowedDataLayouts(_configuration.dataLayoutOptions.value);
  _autoPasContainer->setAllowedNewton3Options(_configuration.newton3Options.value);
  _autoPasContainer->setAllowedTraversals(_configuration.traversalOptions.value);
  _autoPasContainer->setAllowedLoadEstimators(_configuration.loadEstimatorOptions.value);
  _autoPasContainer->setBoxMin(_domainDecomposition.getLocalBoxMin());
  _autoPasContainer->setBoxMax(_domainDecomposition.getLocalBoxMax());
  _autoPasContainer->setCutoff(_configuration.cutoff.value);
  _autoPasContainer->setRelativeOptimumRange(_configuration.relativeOptimumRange.value);
  _autoPasContainer->setMaxTuningPhasesWithoutTest(_configuration.maxTuningPhasesWithoutTest.value);
  _autoPasContainer->setRelativeBlacklistRange(_configuration.relativeBlacklistRange.value);
  _autoPasContainer->setEvidenceFirstPrediction(_configuration.evidenceFirstPrediction.value);
  _autoPasContainer->setExtrapolationMethodOption(_configuration.extrapolationMethodOption.value);
  _autoPasContainer->setNumSamples(_configuration.tuningSamples.value);
  _autoPasContainer->setMaxEvidence(_configuration.tuningMaxEvidence.value);
  _autoPasContainer->setSelectorStrategy(_configuration.selectorStrategy.value);
  _autoPasContainer->setTuningInterval(_configuration.tuningInterval.value);
  _autoPasContainer->setTuningStrategyOption(_configuration.tuningStrategyOption.value);
  _autoPasContainer->setMPIStrategy(_configuration.mpiStrategyOption.value);
  _autoPasContainer->setVerletClusterSize(_configuration.verletClusterSize.value);
  _autoPasContainer->setVerletRebuildFrequency(_configuration.verletRebuildFrequency.value);
  _autoPasContainer->setVerletSkin(_configuration.verletSkinRadius.value);
  _autoPasContainer->setAcquisitionFunction(_configuration.acquisitionFunctionOption.value);
  _autoPasContainer->init();

  // @todo: the object generators should only generate particles relevant for the current rank's domain
  for (auto &particle : _configuration.getParticles()) {
    if (_domainDecomposition.isInsideLocalDomain(particle.getR())) {
      _autoPasContainer->addParticle(particle);
    }
  }

  _configuration.flushParticles();

  if (_configuration.useThermostat.value and _configuration.deltaT.value != 0) {
    if (_configuration.addBrownianMotion.value) {
      Thermostat::addBrownianMotion(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                    _configuration.initTemperature.value);
    }
    Thermostat::apply(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                      _configuration.initTemperature.value, std::numeric_limits<double>::max());
  }

  _timers.initialization.stop();
}

Simulation::~Simulation() { _timers.total.stop(); }

void Simulation::run() {
  const int iterationsPerSuperstep = _configuration.verletRebuildFrequency.value;
  _timers.simulate.start();
  for (int i = 0; i < _configuration.iterations.value; i += iterationsPerSuperstep) {
    executeSupersteps(iterationsPerSuperstep);
  }
  _timers.simulate.stop();

  // Record last state of simulation.
  if (_createVtkFiles) {
    _vtkWriter->recordTimestep(_iteration, *_autoPasContainer, _domainDecomposition);
  }
}

void Simulation::executeSupersteps(const int iterationsPerSuperstep) {
  for (int i = 0; i < iterationsPerSuperstep; ++i) {
    if (_createVtkFiles and _iteration % _configuration.vtkWriteFrequency.value == 0) {
      _timers.vtk.start();
      _vtkWriter->recordTimestep(_iteration, *_autoPasContainer, _domainDecomposition);
      _timers.vtk.stop();
    }

    updatePositions();

    _domainDecomposition.exchangeMigratingParticles(_autoPasContainer);
    _domainDecomposition.exchangeHaloParticles(_autoPasContainer);

    updateForces();

    updateVelocities();

    updateThermostat();

    ++_iteration;

    if (_domainDecomposition.getDomainIndex() == 0) {
      auto [maxIterationsEstimate, maxIterationsIsPrecise] = estimateNumberOfIterations();
      if (not _configuration.dontShowProgressBar.value) {
        printProgress(_iteration, maxIterationsEstimate, maxIterationsIsPrecise);
      }
    }
  }
}

std::tuple<size_t, bool> Simulation::estimateNumberOfIterations() const {
  if (_configuration.tuningPhases.value > 0) {
    // @TODO: this can be improved by considering the tuning strategy
    // This is just a randomly guessed number but seems to fit roughly for default settings.
    size_t configsTestedPerTuningPhase = 90;
    if (_configuration.tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianSearch or
        _configuration.tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianClusterSearch) {
      configsTestedPerTuningPhase = _configuration.tuningMaxEvidence.value;
    }
    auto estimate =
        (_configuration.tuningPhases.value - 1) * _configuration.tuningInterval.value +
        (_configuration.tuningPhases.value * _configuration.tuningSamples.value * configsTestedPerTuningPhase);
    return {estimate, false};
  } else {
    return {_configuration.iterations.value, true};
  }
}

void Simulation::printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise) {
  // percentage of iterations complete
  double fractionDone = static_cast<double>(iterationProgress) / maxIterations;

  // length of the number of maxIterations
  size_t numCharsOfMaxIterations = std::to_string(maxIterations).size();

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
  size_t terminalWidth = getTerminalWidth();

  // the bar should fill the terminal window so subtract everything else (-2 for "] ")
  size_t maxBarWidth = terminalWidth - info.str().size() - progressbar.str().size() - 2;
  // sanity check for underflow
  if (maxBarWidth > terminalWidth) {
    std::cerr << "Warning! Terminal width appears to be too small or could not be read. Disabling progress bar."
              << std::endl;
    _configuration.dontShowProgressBar.value = true;
    return;
  }
  auto barWidth =
      std::max(std::min(static_cast<decltype(maxBarWidth)>(maxBarWidth * (fractionDone)), maxBarWidth), 1ul);
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

std::string Simulation::timerToString(const std::string &name, long timeNS, size_t numberWidth, long maxTime) {
  // only print timers that were actually used
  if (timeNS == 0) {
    return "";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(_floatStringPrecision) << name << " : " << std::setw(numberWidth) << std::right
     << timeNS
     << " ns ("
     // min width of the representation of seconds is numberWidth - 9 (from conversion) + 4 (for dot and digits after)
     << std::setw(numberWidth - 5ul) << ((double)timeNS * 1e-9) << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeNS / (double)maxTime * 100) << "%";
  }
  ss << std::endl;
  return ss.str();
}

void Simulation::updatePositions() {
  _timers.positionUpdate.start();
  TimeDiscretization::calculatePositions(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                         _configuration.deltaT.value);
  _timers.positionUpdate.stop();
}

void Simulation::updateForces() {
  _timers.forceUpdateTotal.start();

  bool isTuningIteration = false;
  _timers.forceUpdatePairwise.start();

  TimeDiscretization::calculatePairwiseForces(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                              _configuration.deltaT.value, _configuration.functorOption.value,
                                              isTuningIteration);

  _timers.forceUpdateTotal.stop();

  auto timeIteration = _timers.forceUpdatePairwise.stop();

  // count time spent for tuning
  if (isTuningIteration) {
    _timers.forceUpdateTuning.addTime(timeIteration);
    ++_numTuningIterations;
  } else {
    _timers.forceUpdateNonTuning.addTime(timeIteration);
    // if the previous iteration was a tuning iteration an the current one is not
    // we have reached the end of a tuning phase
    if (_previousIterationWasTuningIteration) {
      ++_numTuningPhasesCompleted;
    }
  }
  _previousIterationWasTuningIteration = isTuningIteration;

  _timers.forceUpdateTotal.start();
  _timers.forceUpdateGlobal.start();

  if (!_configuration.globalForceIsZero()) {
    TimeDiscretization::calculateGlobalForces(*_autoPasContainer, _configuration.globalForce.value);
  }

  _timers.forceUpdateGlobal.stop();
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

void Simulation::updateThermostat() {
  if (_configuration.useThermostat.value and (_iteration % _configuration.thermostatInterval.value) == 0) {
    _timers.thermostat.start();
    Thermostat::apply(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                      _configuration.targetTemperature.value, _configuration.deltaTemp.value);
    _timers.thermostat.stop();
  }
}
