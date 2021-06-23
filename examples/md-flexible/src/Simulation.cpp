/**
 * @file Simulation.cpp
 * @author J. Körner
 * @date 07.04.2021
 */
#include "Simulation.h"

#include <sys/ioctl.h>
#include <unistd.h>

#include <fstream>
#include <iostream>

#include "BoundaryConditions.h"
#include "Thermostat.h"
#include "TimeDiscretization.h"
#include "configuration/MDFlexConfig.h"
#include "src/ParticleSerializationTools.h"

Simulation::Simulation(const MDFlexConfig &configuration, RegularGridDecomposition &domainDecomposition)
    : _configuration(configuration),
      _domainDecomposition(domainDecomposition),
      _createVtkFiles(not configuration.vtkFileName.value.empty()),
      _vtkWriter(std::make_shared<ParallelVtkWriter>(_configuration.vtkFileName.value, "output", std::to_string(_configuration.iterations.value).size())) {
  _timers.total.start();
  _timers.initialization.start();

  if (_configuration.logFileName.value.empty()) {
    _outputStream = &std::cout;
  } else {
    _logFile = std::make_shared<std::ofstream>();
    _logFile->open(_configuration.logFileName.value);
    _outputStream = &(*_logFile);
  }

  const std::vector<double> localBoxMin = _domainDecomposition.getLocalBoxMin();
  const std::vector<double> localBoxMax = _domainDecomposition.getLocalBoxMax();

  _autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(*_outputStream);
  _autoPasContainer->setAllowedCellSizeFactors(*_configuration.cellSizeFactors.value);
  _autoPasContainer->setAllowedContainers(_configuration.containerOptions.value);
  _autoPasContainer->setAllowedDataLayouts(_configuration.dataLayoutOptions.value);
  _autoPasContainer->setAllowedNewton3Options(_configuration.newton3Options.value);
  _autoPasContainer->setAllowedTraversals(_configuration.traversalOptions.value);
  _autoPasContainer->setAllowedLoadEstimators(_configuration.loadEstimatorOptions.value);
  _autoPasContainer->setBoxMin({localBoxMin[0], localBoxMin[1], localBoxMin[2]});
  _autoPasContainer->setBoxMax({localBoxMax[0], localBoxMax[1], localBoxMax[2]});
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

  // @todo: the object generators should only generate particles relevant for the current ranks domain
  for (auto &particle : _configuration.getParticles()) {
    ParticleType autoPasParticle;
    if (_domainDecomposition.isInsideLocalDomain(particle.position)) {
      autoPasParticle = ParticleSerializationTools::convertParticleAttributesToParticle(particle);
      _autoPasContainer->addParticle(autoPasParticle);
    }
  }

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

void Simulation::run() {
  const int iterationsPerSuperstep = _configuration.verletRebuildFrequency.value;
  for (int i = 0; i < _configuration.iterations.value; i += iterationsPerSuperstep) {
    executeSuperstep(iterationsPerSuperstep);
  }

  // Record last state of simulation.
  if (_createVtkFiles) {
    _vtkWriter->recordTimestep(_iteration, *_autoPasContainer);
  }
}

void Simulation::executeSuperstep(const int iterationsPerSuperstep) {
  for (int i = 0; i < iterationsPerSuperstep; ++i) {

    if(_domainDecomposition.getDomainIndex() == 0){
      auto [maxIterationsEstimate, maxIterationsIsPrecise] = estimateNumberOfIterations();
      if (not _configuration.dontShowProgressBar.value) {
        printProgress(_iteration, maxIterationsEstimate, maxIterationsIsPrecise);
      }
    }

    if (_createVtkFiles and _iteration % _configuration.vtkWriteFrequency.value == 0) {
      _vtkWriter->recordTimestep(_iteration, *_autoPasContainer);
    }

    updatePositions();

    _domainDecomposition.exchangeMigratingParticles(_autoPasContainer);
    _domainDecomposition.exchangeHaloParticles(_autoPasContainer);

    updateForces();

    updateVelocities();

    updateThermostat();

    ++_iteration;
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

bool Simulation::needsMoreIterations() {
  return _iteration < _configuration.iterations.value or _numTuningPhasesCompleted < _configuration.tuningPhases.value;
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
  struct winsize w {};
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  auto terminalWidth = w.ws_col;
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

void Simulation::writeVTKFile() {
  _timers.vtk.start();

  std::string fileBaseName = _configuration.vtkFileName.value;
  // only count number of owned particles here
  const auto numParticles = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned);
  std::ostringstream strstr;
  auto maxNumDigits = std::to_string(_configuration.iterations.value).length();
  strstr << fileBaseName << "_" << getMPISuffix() << std::setfill('0') << std::setw(maxNumDigits) << _iteration
         << ".vtk";
  std::ofstream vtkFile;
  vtkFile.open(strstr.str());

  if (not vtkFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + strstr.str() + "\"");
  }

  vtkFile << "# vtk DataFile Version 2.0\n"
          << "Timestep\n"
          << "ASCII\n";

  // print positions
  vtkFile << "DATASET STRUCTURED_GRID\n"
          << "DIMENSIONS 1 1 1\n"
          << "POINTS " << numParticles << " double\n";

  for (auto iter = _autoPasContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
  }
  vtkFile << "\n";

  vtkFile << "POINT_DATA " << numParticles << "\n";
  // print velocities
  vtkFile << "VECTORS velocities double\n";
  for (auto iter = _autoPasContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto v = iter->getV();
    vtkFile << v[0] << " " << v[1] << " " << v[2] << "\n";
  }
  vtkFile << "\n";

  // print Forces
  vtkFile << "VECTORS forces double\n";
  for (auto iter = _autoPasContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto f = iter->getF();
    vtkFile << f[0] << " " << f[1] << " " << f[2] << "\n";
  }
  vtkFile << "\n";

  // print TypeIDs
  vtkFile << "SCALARS typeIds int\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for (auto iter = _autoPasContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    vtkFile << iter->getTypeId() << "\n";
  }
  vtkFile << "\n";

  // print TypeIDs
  vtkFile << "SCALARS particleIds int\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for (auto iter = _autoPasContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    vtkFile << iter->getID() << "\n";
  }
  vtkFile << "\n";

  vtkFile.close();

  _timers.vtk.stop();
}

std::string Simulation::getMPISuffix() {
  std::string suffix;
#ifdef AUTOPAS_INTERNODE_TUNING
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ostringstream output;
  output << "mpi_rank_" << rank << "_";
  suffix = output.str();
#endif
  return suffix;
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
  TimeDiscretization::calculatePositions(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                         _configuration.deltaT.value);
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
    Thermostat::apply(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                      _configuration.targetTemperature.value, _configuration.deltaTemp.value);
  }
}
