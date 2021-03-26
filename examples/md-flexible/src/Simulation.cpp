/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 01.03.2021
 */

#include "Simulation.h"

#include <sys/ioctl.h>
#include <unistd.h>
#ifdef AUTOPAS_MPI
#include <mpi.h>
#endif

#include <iomanip>
#include <iostream>

#include "BoundaryConditions.h"
#include "Checkpoint.h"
#include "Thermostat.h"
#include "TimeDiscretization.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "autopas/utils/MemoryProfiler.h"

void Simulation::initializeParticlePropertiesLibrary() {
  if (_config->epsilonMap.value.empty()) {
    throw std::runtime_error("No properties found in particle properties library!");
  }

  if (_config->epsilonMap.value.size() != _config->sigmaMap.value.size() or
      _config->epsilonMap.value.size() != _config->massMap.value.size()) {
    throw std::runtime_error("Number of particle properties differ!");
  }

  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibraryType>(_config->cutoff.value);

  for (auto [type, epsilon] : _config->epsilonMap.value) {
    _particlePropertiesLibrary->addType(type, epsilon, _config->sigmaMap.value.at(type),
                                        _config->massMap.value.at(type));
  }
  _particlePropertiesLibrary->calculateMixingCoefficients();
}

void Simulation::initialize(const MDFlexConfig &mdFlexConfig, autopas::AutoPas<ParticleType> &autopas) {
  _timers.init.start();

  _config = std::make_shared<MDFlexConfig>(mdFlexConfig);
  initializeParticlePropertiesLibrary();

  autopas.setAllowedCellSizeFactors(*_config->cellSizeFactors.value);
  autopas.setAllowedContainers(_config->containerOptions.value);
  autopas.setAllowedDataLayouts(_config->dataLayoutOptions.value);
  autopas.setAllowedNewton3Options(_config->newton3Options.value);
  autopas.setAllowedTraversals(_config->traversalOptions.value);
  autopas.setAllowedLoadEstimators(_config->loadEstimatorOptions.value);
  autopas.setBoxMax(_config->boxMax.value);
  autopas.setBoxMin(_config->boxMin.value);
  autopas.setCutoff(_config->cutoff.value);
  autopas.setRelativeOptimumRange(_config->relativeOptimumRange.value);
  autopas.setMaxTuningPhasesWithoutTest(_config->maxTuningPhasesWithoutTest.value);
  autopas.setRelativeBlacklistRange(_config->relativeBlacklistRange.value);
  autopas.setEvidenceFirstPrediction(_config->evidenceFirstPrediction.value);
  autopas.setExtrapolationMethodOption(_config->extrapolationMethodOption.value);
  autopas.setNumSamples(_config->tuningSamples.value);
  autopas.setMaxEvidence(_config->tuningMaxEvidence.value);
  autopas.setSelectorStrategy(_config->selectorStrategy.value);
  autopas.setTuningInterval(_config->tuningInterval.value);
  autopas.setTuningStrategyOption(_config->tuningStrategyOption.value);
  autopas.setMPIStrategy(_config->mpiStrategyOption.value);
  autopas.setVerletClusterSize(_config->verletClusterSize.value);
  autopas.setVerletRebuildFrequency(_config->verletRebuildFrequency.value);
  autopas.setVerletSkin(_config->verletSkinRadius.value);
  autopas.setAcquisitionFunction(_config->acquisitionFunctionOption.value);
  autopas.setOutputSuffix(getMPISuffix());
  autopas::Logger::get()->set_level(_config->logLevel.value);
  autopas.init();

  // load checkpoint
  if (not _config->checkpointfile.value.empty()) {
    Checkpoint::loadParticles(autopas, _config->checkpointfile.value);
  }

  // sanitize simulation end condition
  if (_config->tuningPhases.value > 0) {
    // set iterations to zero because we don't want to consider it
    _config->iterations.value = 0ul;
  }

  // initializing Objects
  for (const auto &object : _config->cubeGridObjects) {
    object.generate(autopas);
  }
  for (const auto &object : _config->cubeGaussObjects) {
    object.generate(autopas);
  }
  for (const auto &object : _config->cubeUniformObjects) {
    object.generate(autopas);
  }
  for (const auto &object : _config->sphereObjects) {
    object.generate(autopas);
  }
  for (const auto &object : _config->cubeClosestPackedObjects) {
    object.generate(autopas);
  }

  // initializing system to initial temperature and Brownian motion
  if (_config->useThermostat.value and _config->deltaT.value != 0) {
    if (_config->addBrownianMotion.value) {
      Thermostat::addBrownianMotion(autopas, *_particlePropertiesLibrary, _config->initTemperature.value);
    }
    // set system to initial temperature
    Thermostat::apply(autopas, *_particlePropertiesLibrary, _config->initTemperature.value,
                      std::numeric_limits<double>::max());
  }

  _timers.init.stop();
}

template <class FunctorType>
void Simulation::calculateForces(autopas::AutoPas<ParticleType> &autopas) {
  _timers.forceUpdateTotal.start();

  // pairwise forces

  _timers.forceUpdatePairwise.start();

  FunctorType functor{autopas.getCutoff(), *_particlePropertiesLibrary};
  bool tuningIteration = autopas.iteratePairwise(&functor);

  _timers.forceUpdateTotal.stop();
  auto timeIteration = _timers.forceUpdatePairwise.stop();

  // count time spent for tuning
  if (tuningIteration) {
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
  _previousIterationWasTuningIteration = tuningIteration;

  // global forces

  _timers.forceUpdateTotal.start();
  _timers.forceUpdateGlobal.start();
  globalForces(autopas);
  _timers.forceUpdateGlobal.stop();
  _timers.forceUpdateTotal.stop();
}

void Simulation::globalForces(autopas::AutoPas<ParticleType> &autopas) {
  // skip application of zero force
  if (_config->globalForceIsZero()) {
    return;
  }

#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(autopas)
#endif
  for (auto particleItr = autopas.begin(autopas::IteratorBehavior::ownedOnly); particleItr.isValid(); ++particleItr) {
    particleItr->addF(_config->globalForce.value);
  }
}

void Simulation::simulate(autopas::AutoPas<ParticleType> &autopas) {
  this->_homogeneity = Simulation::calculateHomogeneity(autopas);
  _timers.simulate.start();

  auto [maxIterationsEstimate, maxIterationsIsPrecise] = estimateNumIterations();

  // main simulation loop
  for (; needsMoreIterations(); ++_iteration) {
    if (not _config->dontShowProgressBar.value) {
      printProgress(_iteration, maxIterationsEstimate, maxIterationsIsPrecise);
    }

    std::cout << "This is " << getMPISuffix() << ". Homogeneity: " << calculateHomogeneity(autopas) << "\n";

    // only do time step related stuff when there actually is time-stepping
    if (_config->deltaT.value != 0) {
      // only write vtk files periodically and if a filename is given.
      if ((not _config->vtkFileName.value.empty()) and _iteration % _config->vtkWriteFrequency.value == 0) {
        this->writeVTKFile(autopas);
      }

      // calculate new positions
      _timers.positionUpdate.start();
      TimeDiscretization::calculatePositions(autopas, *_particlePropertiesLibrary, _config->deltaT.value);
      _timers.positionUpdate.stop();

      // apply boundary conditions AFTER the position update!
      if (_config->periodic.value) {
        _timers.boundaries.start();
        BoundaryConditions::applyPeriodic(autopas, false);
        _timers.boundaries.stop();
      } else {
        throw std::runtime_error(
            "Simulation::simulate(): at least one boundary condition has to be set. Please enable the periodic "
            "boundary conditions!");
      }
    }
    // invoke the force calculation with the functor specified in the configuration
    switch (this->_config->functorOption.value) {
      case MDFlexConfig::FunctorOption::lj12_6: {
        this->calculateForces<autopas::LJFunctor<ParticleType, _shifting, _mixing>>(autopas);
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_Globals: {
        this->calculateForces<
            autopas::LJFunctor<ParticleType, _shifting, _mixing, autopas::FunctorN3Modes::Both, /* globals */ true>>(
            autopas);
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_AVX: {
        this->calculateForces<autopas::LJFunctorAVX<ParticleType, _shifting, _mixing>>(autopas);
        break;
      }
    }
    // only show memory usage in when the logger is set to debug
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << std::endl;
    }

    // only do time step related stuff when there actually is time-stepping
    if (_config->deltaT.value != 0) {
      _timers.velocityUpdate.start();
      TimeDiscretization::calculateVelocities(autopas, *_particlePropertiesLibrary, _config->deltaT.value);
      _timers.velocityUpdate.stop();

      // applying Velocity scaling with Thermostat:
      if (_config->useThermostat.value and (_iteration % _config->thermostatInterval.value) == 0) {
        _timers.thermostat.start();
        Thermostat::apply(autopas, *_particlePropertiesLibrary, _config->targetTemperature.value,
                          _config->deltaTemp.value);
        _timers.thermostat.stop();
      }
    }
  }
  // final update for a full progress bar
  if (not _config->dontShowProgressBar.value) {
    // The last update is precise, so we know the number of iterations.
    printProgress(_iteration, _iteration, true);
    // The progress bar does not end the line. Since this is the last progress bar, end the line here.
    std::cout << std::endl;
  }

  // update temperature for generated config output
  if (_config->useThermostat.value) {
    _timers.thermostat.start();
    _config->initTemperature.value = Thermostat::calcTemperature(autopas, *_particlePropertiesLibrary);
    _timers.thermostat.stop();
  }

  // writes final state of the simulation
  if ((not _config->vtkFileName.value.empty())) {
    _timers.boundaries.start();
    BoundaryConditions::applyPeriodic(autopas, true);
    _timers.boundaries.stop();
    this->writeVTKFile(autopas);
  }

  _timers.simulate.stop();
}

void Simulation::printStatistics(autopas::AutoPas<ParticleType> &autopas) {
  using namespace std;
  size_t flopsPerKernelCall;

  switch (_config->functorOption.value) {
    case MDFlexConfig::FunctorOption ::lj12_6: {
      flopsPerKernelCall = autopas::LJFunctor<ParticleType, _shifting, _mixing>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_Globals: {
      flopsPerKernelCall = autopas::LJFunctor<ParticleType, _shifting, _mixing, autopas::FunctorN3Modes::Both,
                                              /* globals */ true>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_AVX: {
      flopsPerKernelCall = autopas::LJFunctorAVX<ParticleType, _shifting, _mixing>::getNumFlopsPerKernelCall();
      break;
    }
    default:
      throw std::runtime_error("Invalid Functor choice");
  }

  auto durationTotal = _timers.total.stop();
  auto durationSimulate = _timers.simulate.getTotalTime();
  auto durationSimulateSec = durationSimulate * 1e-9;

  // take total time as base for formatting since this should be the longest
  auto digitsTimeTotalNS = std::to_string(durationTotal).length();

  // Statistics
  cout << endl;
  cout << "Total number of particles at end of Simulation: "
       << autopas.getNumberOfParticles(autopas::IteratorBehavior::haloAndOwned) << endl;
  cout << "  Owned: " << autopas.getNumberOfParticles(autopas::IteratorBehavior::ownedOnly) << endl;
  cout << "  Halo : " << autopas.getNumberOfParticles(autopas::IteratorBehavior::haloOnly) << endl;
  cout << "Standard Deviation of Homogeneity    : " << _homogeneity << endl;

  cout << fixed << setprecision(_floatStringPrecision);
  cout << "Measurements:" << endl;
  cout << timerToString("Time total      ", durationTotal, digitsTimeTotalNS);
  cout << timerToString("  Initialization", _timers.init.getTotalTime(), digitsTimeTotalNS, durationTotal);
  cout << timerToString("  Simulation    ", durationSimulate, digitsTimeTotalNS, durationTotal);
  cout << timerToString("    Boundaries  ", _timers.boundaries.getTotalTime(), digitsTimeTotalNS, durationSimulate);
  cout << timerToString("    Position    ", _timers.positionUpdate.getTotalTime(), digitsTimeTotalNS, durationSimulate);
  cout << timerToString("    Force       ", _timers.forceUpdateTotal.getTotalTime(), digitsTimeTotalNS,
                        durationSimulate);
  cout << timerToString("      Tuning    ", _timers.forceUpdateTuning.getTotalTime(), digitsTimeTotalNS,
                        _timers.forceUpdateTotal.getTotalTime());
  cout << timerToString("      NonTuning ", _timers.forceUpdateNonTuning.getTotalTime(), digitsTimeTotalNS,
                        _timers.forceUpdateTotal.getTotalTime());
  cout << timerToString("    Velocity    ", _timers.velocityUpdate.getTotalTime(), digitsTimeTotalNS, durationSimulate);
  cout << timerToString("    VTK         ", _timers.vtk.getTotalTime(), digitsTimeTotalNS, durationSimulate);
  cout << timerToString("    Thermostat  ", _timers.thermostat.getTotalTime(), digitsTimeTotalNS, durationSimulate);

  cout << timerToString("One iteration   ", _timers.simulate.getTotalTime() / _iteration, digitsTimeTotalNS,
                        durationTotal);
  auto mfups = autopas.getNumberOfParticles(autopas::IteratorBehavior::ownedOnly) * _iteration * 1e-6 /
               (_timers.forceUpdateTotal.getTotalTime() * 1e-9);  // 1e-9 for ns to s, 1e-6 for M in MFUP
  cout << "Tuning iterations: " << _numTuningIterations << " / " << _iteration << " = "
       << ((double)_numTuningIterations / _iteration * 100) << "%" << endl;
  cout << "MFUPs/sec    : " << mfups << endl;

  if (_config->dontMeasureFlops.value) {
    autopas::FlopCounterFunctor<ParticleType> flopCounterFunctor(autopas.getCutoff());
    autopas.iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * _iteration;
    // approximation for flops of verlet list generation
    if (autopas.getContainerType() == autopas::ContainerOption::verletLists)
      flops += flopCounterFunctor.getDistanceCalculations() *
               decltype(flopCounterFunctor)::numFlopsPerDistanceCalculation *
               floor(_iteration / _config->verletRebuildFrequency.value);

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationSimulateSec << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }
}

const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &Simulation::getPpl() const {
  return _particlePropertiesLibrary;
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

bool Simulation::needsMoreIterations() const {
  return _iteration < _config->iterations.value or _numTuningPhasesCompleted < _config->tuningPhases.value;
}

std::string Simulation::getMPISuffix() const {
  std::string suffix;
#ifdef AUTOPAS_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ostringstream output;
  output << "mpi_rank_" << rank << "_";
  suffix = output.str();
#endif
  return suffix;
}

void Simulation::writeVTKFile(autopas::AutoPas<ParticleType> &autopas) {
  _timers.vtk.start();

  std::string fileBaseName = _config->vtkFileName.value;
  // only count number of owned particles here
  const auto numParticles = autopas.getNumberOfParticles(autopas::IteratorBehavior::ownedOnly);
  std::ostringstream strstr;
  auto maxNumDigits = std::to_string(_config->iterations.value).length();
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

  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
  }
  vtkFile << "\n";

  vtkFile << "POINT_DATA " << numParticles << "\n";
  // print velocities
  vtkFile << "VECTORS velocities double\n";
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto v = iter->getV();
    vtkFile << v[0] << " " << v[1] << " " << v[2] << "\n";
  }
  vtkFile << "\n";

  // print Forces
  vtkFile << "VECTORS forces double\n";
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto f = iter->getF();
    vtkFile << f[0] << " " << f[1] << " " << f[2] << "\n";
  }
  vtkFile << "\n";

  // print TypeIDs
  vtkFile << "SCALARS typeIds int\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    vtkFile << iter->getTypeId() << "\n";
  }
  vtkFile << "\n";

  // print TypeIDs
  vtkFile << "SCALARS particleIds int\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    vtkFile << iter->getID() << "\n";
  }
  vtkFile << "\n";

  vtkFile.close();

  _timers.vtk.stop();
}

double Simulation::calculateHomogeneity(autopas::AutoPas<ParticleType> &autopas) const {
  size_t numberOfParticles = autopas.getNumberOfParticles();
  // approximately the resolution we want to get.
  size_t numberOfCells = ceil(numberOfParticles / 10.);

  std::array<double, 3> startCorner = autopas.getBoxMin();
  std::array<double, 3> endCorner = autopas.getBoxMax();
  std::array<double, 3> domainSizePerDimension = {};
  for (int i = 0; i < 3; ++i) {
    domainSizePerDimension[i] = endCorner[i] - startCorner[i];
  }

  // get cellLength which is equal in each direction, derived from the domainsize and the requested number of cells
  double volume = domainSizePerDimension[0] * domainSizePerDimension[1] * domainSizePerDimension[2];
  double cellVolume = volume / numberOfCells;
  double cellLength = cbrt(cellVolume);

  // calculate the size of the boundary cells, which might be smaller then the other cells
  std::array<size_t, 3> cellsPerDimension = {};
  // size of the last cell layer per dimension. This cell might get truncated to fit in the domain.
  std::array<double, 3> outerCellSizePerDimension = {};
  for (int i = 0; i < 3; ++i) {
    outerCellSizePerDimension[i] =
        domainSizePerDimension[i] - (floor(domainSizePerDimension[i] / cellLength) * cellLength);
    cellsPerDimension[i] = ceil(domainSizePerDimension[i] / cellLength);
  }
  // Actual number of cells we end up with
  numberOfCells = cellsPerDimension[0] * cellsPerDimension[1] * cellsPerDimension[2];

  std::vector<size_t> particlesPerCell(numberOfCells, 0);
  std::vector<double> allVolumes(numberOfCells, 0);

  // add particles accordingly to their cell to get the amount of particles in each cell
  for (auto particleItr = autopas.begin(autopas::IteratorBehavior::ownedOnly); particleItr.isValid(); ++particleItr) {
    std::array<double, 3> particleLocation = particleItr->getR();
    std::array<size_t, 3> index = {};
    for (int i = 0; i < particleLocation.size(); i++) {
      index[i] = particleLocation[i] / cellLength;
    }
    const size_t cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(index, cellsPerDimension);
    particlesPerCell[cellIndex] += 1;
    // calculate the size of the current cell
    allVolumes[cellIndex] = 1;
    for (int i = 0; i < cellsPerDimension.size(); ++i) {
      // the last cell layer has a special size
      if (index[i] == cellsPerDimension[i] - 1) {
        allVolumes[cellIndex] *= outerCellSizePerDimension[i];
      } else {
        allVolumes[cellIndex] *= cellLength;
      }
    }
  }

  // calculate density for each cell
  std::vector<double> densityPerCell(numberOfCells, 0.0);
  for (int i = 0; i < particlesPerCell.size(); i++) {
    densityPerCell[i] =
        (allVolumes[i] == 0) ? 0 : (particlesPerCell[i] / allVolumes[i]);  // make sure there is no division of zero
  }

  // get mean and reserve variable for variance
  double mean = numberOfParticles / volume;
  double variance = 0.0;

  // calculate variance
  for (int r = 0; r < densityPerCell.size(); ++r) {
    double distance = densityPerCell[r] - mean;
    variance += (distance * distance / densityPerCell.size());
  }

  // finally calculate standard deviation
  return sqrt(variance);
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
    _config->dontShowProgressBar.value = true;
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

std::tuple<size_t, bool> Simulation::estimateNumIterations() const {
  if (_config->tuningPhases.value > 0) {
    // @TODO: this can be improved by considering the tuning strategy
    // This is just a randomly guessed number but seems to fit roughly for default settings.
    size_t configsTestedPerTuningPhase = 90;
    if (_config->tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianSearch or
        _config->tuningStrategyOption.value == autopas::TuningStrategyOption::bayesianClusterSearch) {
      configsTestedPerTuningPhase = _config->tuningMaxEvidence.value;
    }
    auto estimate = (_config->tuningPhases.value - 1) * _config->tuningInterval.value +
                    (_config->tuningPhases.value * _config->tuningSamples.value * configsTestedPerTuningPhase);
    return {estimate, false};
  } else {
    return {_config->iterations.value, true};
  }
}
