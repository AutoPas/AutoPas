/**
 * @file Simulation.h
 * @author N. Fottner
 * @date 12/05/19
 */
#pragma once

#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "BoundaryConditions.h"
#include "Checkpoint.h"
#include "Generator.h"
#include "PrintableMolecule.h"
#include "Thermostat.h"
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/utils/MemoryProfiler.h"
#include "parsing/MDFlexConfig.h"

template <class Particle, class ParticleCell>
class Simulation {
 public:
  /**
   * Constructor.
   *
   * This starts the timer for total simulation time.
   */
  explicit Simulation() { _timers.total.start(); };

  /**
   * Destructor.
   *
   * Closes the log file if applicable.
   */
  ~Simulation() {
    if (not _config->logFileName.empty()) {
      _logFile.close();
    }
  }

  /**
   * Writes a VTK file for the current state of the AutoPas object
   * @tparam AutoPasTemplate Template for the templetized autopas type.
   * @param filename
   * @param numParticles
   * @param autopas
   */
  void writeVTKFile(unsigned int iteration) {
    _timers.vtk.start();

    std::string fileBaseName = _config->vtkFileName;
    // as _autopas.getNumberOfParticles return number of haloAndOwned Particles, we need number of owned Particles
    const auto numParticles = this->getNumParticles();
    std::ostringstream strstr;
    auto maxNumDigits = std::to_string(_config->iterations).length();
    strstr << fileBaseName << "_" << std::setfill('0') << std::setw(maxNumDigits) << iteration << ".vtk";
    std::ofstream vtkFile;
    vtkFile.open(strstr.str());

    vtkFile << "# vtk DataFile Version 2.0" << std::endl;
    vtkFile << "Timestep" << std::endl;
    vtkFile << "ASCII" << std::endl;

    // print positions
    vtkFile << "DATASET STRUCTURED_GRID" << std::endl;
    vtkFile << "DIMENSIONS 1 1 1" << std::endl;
    vtkFile << "POINTS " << numParticles << " double" << std::endl;
    for (auto iter = _autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      auto pos = iter->getR();
      vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    }
    vtkFile << std::endl;

    vtkFile << "POINT_DATA " << numParticles << std::endl;
    // print velocities
    vtkFile << "VECTORS velocities double" << std::endl;
    for (auto iter = _autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      auto v = iter->getV();
      vtkFile << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }
    vtkFile << std::endl;

    // print Forces
    vtkFile << "VECTORS forces double" << std::endl;
    for (auto iter = _autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      auto f = iter->getF();
      vtkFile << f[0] << " " << f[1] << " " << f[2] << std::endl;
    }
    vtkFile << std::endl;

    vtkFile.close();

    _timers.vtk.stop();
  }

  /**
   * Initialized the ParticlePropertiesLibrary for usage in functor
   * during Simualtion::initialize call
   * with the parsed values in the yamlParser
   */
  void initializeParticlePropertiesLibrary();

  /**
   * Initializes the AutoPas Object with the given config and initializes the simulation domain with the Object
   * Generators.
   */
  void initialize(const MDFlexConfig &mdFlexConfig);

  /**
   * Does the ForceCalculation
   * @tparam Force Calculation Functor
   */
  template <class FunctorType>
  void calculateForces();

  /**
   * This function processes the main simulation loop
   * -calls the time discretization class(calculate fores, etc ...)
   * -do the output each timestep
   * -collects the duration of every Calculation(Position,Force,Velocity)
   */
  void simulate();

  /**
   * Getter for AutoPas object.
   * @return Autopas Object
   */
  autopas::AutoPas<Particle, ParticleCell> *getAutopas() const;

  /**
   * Return the current number of owned Particles in the AutoPas Object.
   * @return
   */
  size_t getNumParticles() {
    size_t numberOfParticles = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
    for (auto iter = _autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      numberOfParticles++;
    }
    return numberOfParticles;
  }

  /**
   * Prints statistics like duration of calculation etc of the Simulation.
   */
  void printStatistics();

  /**
   * Getter for ParticlePropertiesLibrary of Simulation
   * @return unique_prt(ParticlePropertiesLibrary)
   */
  const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &getPpl() const;

 private:
  autopas::AutoPas<Particle, ParticleCell> _autopas;
  std::shared_ptr<MDFlexConfig> _config;
  std::ofstream _logFile;
  std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> _particlePropertiesLibrary;
  std::unique_ptr<Thermostat<decltype(_autopas), std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>>
      _thermostat;
  std::unique_ptr<
      TimeDiscretization<decltype(_autopas), std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>>
      _timeDiscretization;

  struct timers {
    autopas::utils::Timer positionUpdate, forceUpdate, velocityUpdate, simulate, vtk, init, total, thermostat,
        boundaries;
  } _timers;

  /**
   * Precision of floating point numbers printed.
   */
  constexpr static auto _floatPrecision = 3;

  /**
   * Convert a time and a name to a propperly formatted string.
   * @param name incl. offset.
   * @param timeMS in microseconds.
   * @param numberWidth Width to which the time should be offset.
   * @param maxTime if passed the percentage of timeMS of maxTime is appended.
   * @return formatted std::string
   */
  std::string timerToString(std::string name, long timeMS, size_t numberWidth = 0, long maxTime = 0);
};

template <class Particle, class ParticleCell>
autopas::AutoPas<Particle, ParticleCell> *Simulation<Particle, ParticleCell>::getAutopas() const {
  return _autopas.get();
}

template <typename Particle, typename ParticleCell>
void Simulation<Particle, ParticleCell>::initializeParticlePropertiesLibrary() {
  if (_config->epsilonMap.empty()) {
    throw std::runtime_error("No properties found in particle properties library!");
  }

  if (_config->epsilonMap.size() != _config->sigmaMap.size() or _config->epsilonMap.size() != _config->massMap.size()) {
    throw std::runtime_error("Number of particle properties differ!");
  }

  _particlePropertiesLibrary = std::make_unique<std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>();

  for (auto [type, epsilon] : _config->epsilonMap) {
    _particlePropertiesLibrary->addType(type, epsilon, _config->sigmaMap.at(type), _config->massMap.at(type));
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(const MDFlexConfig &mdFlexConfig) {
  _timers.init.start();

  _config = std::make_shared<MDFlexConfig>(mdFlexConfig);
  initializeParticlePropertiesLibrary();
  _timeDiscretization = std::make_unique<
      TimeDiscretization<decltype(_autopas), std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>>(
      _config->deltaT, *_particlePropertiesLibrary);
  auto logFileName(_config->logFileName);
  auto verletRebuildFrequency(_config->verletRebuildFrequency);
  auto logLevel(_config->logLevel);
  auto &cellSizeFactors(_config->cellSizeFactors);
  auto tuningStrategy(_config->tuningStrategyOption);
  auto containerChoice(_config->containerOptions);
  auto selectorStrategy(_config->selectorStrategy);
  auto cutoff(_config->cutoff);
  auto dataLayoutOptions(_config->dataLayoutOptions);
  auto newton3Options(_config->newton3Options);
  auto traversalOptions(_config->traversalOptions);
  auto tuningInterval(_config->tuningInterval);
  auto tuningSamples(_config->tuningSamples);
  auto verletSkinRadius(_config->verletSkinRadius);
  auto cubesGrid(_config->cubeGridObjects);
  auto cubesGauss(_config->cubeGaussObjects);
  auto cubesUniform(_config->cubeUniformObjects);
  auto spheres(_config->sphereObjects);

  // select either std::out or a logfile for autopas log output.
  // This does not affect md-flex output.
  std::streambuf *streamBuf;
  if (logFileName.empty()) {
    streamBuf = std::cout.rdbuf();
  } else {
    _logFile.open(logFileName);
    streamBuf = _logFile.rdbuf();
  }
  std::ostream outputStream(streamBuf);
  _autopas.setAllowedCellSizeFactors(*cellSizeFactors);
  _autopas.setAllowedContainers(containerChoice);
  _autopas.setAllowedDataLayouts(dataLayoutOptions);
  _autopas.setAllowedNewton3Options(newton3Options);
  _autopas.setAllowedTraversals(traversalOptions);
  _autopas.setBoxMax(_config->boxMax);
  _autopas.setBoxMin(_config->boxMin);
  _autopas.setCutoff(cutoff);
  _autopas.setNumSamples(tuningSamples);
  _autopas.setSelectorStrategy(selectorStrategy);
  _autopas.setTuningInterval(tuningInterval);
  _autopas.setTuningStrategyOption(tuningStrategy);
  _autopas.setVerletClusterSize(_config->verletClusterSize);
  _autopas.setVerletRebuildFrequency(_config->verletRebuildFrequency);
  _autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  _autopas.setVerletSkin(verletSkinRadius);
  autopas::Logger::get()->set_level(logLevel);
  _autopas.init();

  // load checkpoint
  if (not _config->checkpointfile.empty()) {
    Checkpoint<decltype(_autopas)>::loadParticles(_autopas, _config->checkpointfile);
  }

  // initializing Objects
  for (const auto &grid : cubesGrid) {
    Generator::cubeGrid<Particle, ParticleCell>(_autopas, grid);
  }
  for (const auto &cube : cubesGauss) {
    Generator::cubeGauss<Particle, ParticleCell>(_autopas, cube);
  }
  for (const auto &cube : cubesUniform) {
    Generator::cubeRandom<Particle, ParticleCell>(_autopas, cube);
  }
  for (const auto &sphere : spheres) {
    Generator::sphere<Particle, ParticleCell>(_autopas, sphere);
  }

  // initilizing Thermostat
  if (_config->useThermostat) {
    _thermostat = std::make_unique<
        Thermostat<decltype(_autopas), std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>>(
        _config->initTemperature, _config->targetTemperature, _config->deltaTemp, *_particlePropertiesLibrary);
    _thermostat->addBrownianMotion(_autopas, _config->useCurrentTempForBrownianMotion);
  }

  _timers.init.stop();
}

template <class Particle, class ParticleCell>
template <class FunctorType>
void Simulation<Particle, ParticleCell>::calculateForces() {
  _timers.forceUpdate.start();

  auto functor = FunctorType(_autopas.getCutoff(), 0.0, *_particlePropertiesLibrary);
  _autopas.iteratePairwise(&functor);

  _timers.forceUpdate.stop();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate() {
  _timers.simulate.start();

  // main simulation loop
  for (size_t iteration = 0; iteration < _config->iterations; ++iteration) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Iteration " << iteration << std::endl;
    }

    // only write vtk files periodically and if a filename is given
    if ((not _config->vtkFileName.empty()) and iteration % _config->vtkWriteFrequency == 0) {
      this->writeVTKFile(iteration);
    }

    if (_config->periodic) {
      _timers.boundaries.start();
      BoundaryConditions<ParticleCell>::applyPeriodic(_autopas);
      _timers.boundaries.stop();
    }
    _timers.positionUpdate.start();
    _timeDiscretization->calculatePositions(_autopas);
    _timers.positionUpdate.stop();

    switch (this->_config->functorOption) {
      case MDFlexConfig::FunctorOption::lj12_6: {
        this->calculateForces<autopas::LJFunctor<Particle, ParticleCell, /* mixing */ true>>();
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_Globals: {
        this->calculateForces<autopas::LJFunctor<Particle, ParticleCell, /* mixing */ true,
                                                 autopas::FunctorN3Modes::Both, /* globals */ true>>();
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_AVX: {
        this->calculateForces<autopas::LJFunctorAVX<Particle, ParticleCell, /* mixing */ true>>();
        break;
      }
    }
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << std::endl;
    }
    _timers.velocityUpdate.start();
    _timeDiscretization->calculateVelocities(_autopas);
    _timers.velocityUpdate.stop();

    // applying Velocity scaling with Thermostat:
    if (_config->useThermostat and (iteration % _config->thermostatInterval) == 0) {
      _timers.thermostat.start();
      _thermostat->apply(_autopas);
      _timers.thermostat.stop();
    }
  }

  // writes final state of the simulation
  if ((not _config->vtkFileName.empty())) {
    this->writeVTKFile(_config->iterations);
  }

  _timers.simulate.stop();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::printStatistics() {
  using namespace std;
  size_t flopsPerKernelCall;

  switch (_config->functorOption) {
    case MDFlexConfig::FunctorOption ::lj12_6: {
      flopsPerKernelCall = autopas::LJFunctor<Particle, ParticleCell, /* mixing */ true>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_Globals: {
      flopsPerKernelCall = autopas::LJFunctor<Particle, ParticleCell, /* mixing */ true, autopas::FunctorN3Modes::Both,
                                              /* globals */ true>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_AVX: {
      flopsPerKernelCall = autopas::LJFunctorAVX<Particle, ParticleCell, /* mixing */ true>::getNumFlopsPerKernelCall();
      break;
    }
    default:
      throw std::runtime_error("Invalid Functor choice");
  }

  auto durationTotal = _timers.total.stop();
  auto durationSimulate = _timers.simulate.getTotalTime();
  auto durationSimulateSec = durationSimulate * 1e-6;

  // take total time as base for formatting since this should be the longest
  auto digitsTimeTotalMuS = std::to_string(durationTotal).length();

  // time statistics
  //  cout << "Simulation duration without initilization: " << _timers.durationSimulate << " \u03bcs" << endl;
  // Statistics
  cout << fixed << setprecision(_floatPrecision);
  cout << endl << "Measurements:" << endl;
  cout << timerToString("Time total      ", durationTotal, digitsTimeTotalMuS) << endl;
  cout << timerToString("  Initialization", _timers.init.getTotalTime(), digitsTimeTotalMuS, durationTotal) << endl;
  cout << timerToString("  Simulation    ", durationSimulate, digitsTimeTotalMuS, durationTotal) << endl;
  if (_config->periodic) {
    cout << timerToString("    Boundaries  ", _timers.boundaries.getTotalTime(), digitsTimeTotalMuS, durationSimulate)
         << endl;
  }
  cout << timerToString("    Position    ", _timers.positionUpdate.getTotalTime(), digitsTimeTotalMuS, durationSimulate)
       << endl;
  cout << timerToString("    Force       ", _timers.forceUpdate.getTotalTime(), digitsTimeTotalMuS, durationSimulate)
       << endl;
  cout << timerToString("    Velocity    ", _timers.velocityUpdate.getTotalTime(), digitsTimeTotalMuS, durationSimulate)
       << endl;
  if (not _config->vtkFileName.empty()) {
    cout << timerToString("    VTK         ", _timers.vtk.getTotalTime(), digitsTimeTotalMuS, durationSimulate) << endl;
  }
  if (_config->useThermostat) {
    cout << timerToString("    Thermostat  ", _timers.thermostat.getTotalTime(), digitsTimeTotalMuS, durationSimulate)
         << endl;
  }

  auto numIterations = _config->iterations;

  if (numIterations > 0) {
    cout << timerToString("One iteration   ", _timers.simulate.getTotalTime() / numIterations, digitsTimeTotalMuS,
                          durationTotal)
         << endl;
  }
  auto mfups = _autopas.getNumberOfParticles() * numIterations / durationSimulateSec * 1e-6;
  cout << "MFUPs/sec    : " << mfups << endl;

  if (_config->measureFlops) {
    autopas::FlopCounterFunctor<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> flopCounterFunctor(
        _autopas.getCutoff());
    _autopas.iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * numIterations;
    // approximation for flops of verlet list generation
    if (_autopas.getContainerType() == autopas::ContainerOption::verletLists)
      flops +=
          flopCounterFunctor.getDistanceCalculations() *
          autopas::FlopCounterFunctor<PrintableMolecule,
                                      autopas::FullParticleCell<PrintableMolecule>>::numFlopsPerDistanceCalculation *
          floor(numIterations / _config->verletRebuildFrequency);

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationSimulateSec << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }
}

template <class Particle, class ParticleCell>
const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &Simulation<Particle, ParticleCell>::getPpl() const {
  return _particlePropertiesLibrary;
}
template <class Particle, class ParticleCell>

std::string Simulation<Particle, ParticleCell>::timerToString(std::string name, long timeMS, size_t numberWidth,
                                                              long maxTime) {
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(_floatPrecision);
  ss << name << " : " << std::setw(numberWidth) << std::right << timeMS << " \u03bcs (" << ((double)timeMS * 1e-6)
     << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeMS / (double)maxTime * 100) << "%";
  }
  return ss.str();
}
