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
  ~Simulation() = default;

  /**
   * Writes a VTK file for the current state of the AutoPas object.
   * @tparam AutoPasTemplate Template for the templetized autopas type.
   * @param filename
   * @param numParticles
   * @param autopas
   */
  void writeVTKFile(unsigned int iteration) {
    _timers.vtk.start();

    std::string fileBaseName = _config->vtkFileName;
    // only count number of owned particles here
    const auto numParticles = this->_autopas.getNumberOfParticles(autopas::IteratorBehavior::ownedOnly);
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

    // print TypeIDs
    vtkFile << "SCALARS typeIds int" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (auto iter = _autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      vtkFile << iter->getTypeId() << std::endl;
    }
    vtkFile << std::endl;

    vtkFile.close();

    _timers.vtk.stop();
  }

  /**
   * Initializes the ParticlePropertiesLibrary with properties from _config.
   */
  void initializeParticlePropertiesLibrary();

  /**
   * Initializes the AutoPas Object with the given config and initializes the simulation domain with the Object
   * Generators.
   */
  void initialize(const MDFlexConfig &mdFlexConfig);

  /**
   * Calculates the pairwise forces in the system and measures the runtime.
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
   * Prints statistics like duration of calculation etc of the Simulation.
   */
  void printStatistics();

  /**
   * Getter for ParticlePropertiesLibrary of Simulation.
   * @return unique_prt(ParticlePropertiesLibrary)
   */
  const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &getPpl() const;

 private:
  using AutoPasType = autopas::AutoPas<Particle, ParticleCell>;
  using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<double, size_t>;

  AutoPasType _autopas;
  std::shared_ptr<MDFlexConfig> _config;
  std::ofstream _logFile;
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;

  struct timers {
    autopas::utils::Timer positionUpdate, forceUpdateTotal, forceUpdateTuning, forceUpdateNonTuning, velocityUpdate,
        simulate, vtk, init, total, thermostat, boundaries;
  } _timers;

  /**
   * Precision of floating point numbers printed.
   */
  constexpr static auto _floatStringPrecision = 3;

  /**
   * Convert a time and a name to a properly formatted string.
   * @param name incl. offset.
   * @param timeMS in microseconds.
   * @param numberWidth Width to which the time should be offset.
   * @param maxTime if passed the percentage of timeMS of maxTime is appended.
   * @return formatted std::string
   */
  std::string timerToString(const std::string &name, long timeMS, size_t numberWidth = 0, long maxTime = 0);
};

template <typename Particle, typename ParticleCell>
void Simulation<Particle, ParticleCell>::initializeParticlePropertiesLibrary() {
  if (_config->epsilonMap.empty()) {
    throw std::runtime_error("No properties found in particle properties library!");
  }

  if (_config->epsilonMap.size() != _config->sigmaMap.size() or _config->epsilonMap.size() != _config->massMap.size()) {
    throw std::runtime_error("Number of particle properties differ!");
  }

  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibraryType>();

  for (auto [type, epsilon] : _config->epsilonMap) {
    _particlePropertiesLibrary->addType(type, epsilon, _config->sigmaMap.at(type), _config->massMap.at(type));
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(const MDFlexConfig &mdFlexConfig) {
  _timers.init.start();

  _config = std::make_shared<MDFlexConfig>(mdFlexConfig);
  initializeParticlePropertiesLibrary();
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
    Checkpoint::loadParticles(_autopas, _config->checkpointfile);
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

  // initializing system to initial temperature and Brownian motion
  if (_config->useThermostat and _config->deltaT != 0) {
    if (_config->addBrownianMotion) {
      Thermostat::addBrownianMotion(_autopas, *_particlePropertiesLibrary, _config->initTemperature);
    }
    // set system to initial temperature
    Thermostat::apply(_autopas, *_particlePropertiesLibrary, _config->initTemperature,
                      std::numeric_limits<double>::max());
  }

  _timers.init.stop();
}

template <class Particle, class ParticleCell>
template <class FunctorType>
void Simulation<Particle, ParticleCell>::calculateForces() {
  _timers.forceUpdateTotal.start();

  FunctorType functor{_autopas.getCutoff(), 0.0, *_particlePropertiesLibrary};
  bool tuningIteration = _autopas.iteratePairwise(&functor);

  auto timeIteration = _timers.forceUpdateTotal.stop();
  if (tuningIteration) {
    _timers.forceUpdateTuning.addTime(timeIteration);
  } else {
    _timers.forceUpdateNonTuning.addTime(timeIteration);
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate() {
  _timers.simulate.start();

  // main simulation loop
  for (size_t iteration = 0; iteration < _config->iterations; ++iteration) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Iteration " << iteration << std::endl;
    }

    if (_config->deltaT != 0) {
      // only write vtk files periodically and if a filename is given.
      if ((not _config->vtkFileName.empty()) and iteration % _config->vtkWriteFrequency == 0) {
        this->writeVTKFile(iteration);
      }

      if (_config->periodic) {
        _timers.boundaries.start();
        BoundaryConditions<ParticleCell>::applyPeriodic(_autopas);
        _timers.boundaries.stop();
      }
      _timers.positionUpdate.start();
      TimeDiscretization::calculatePositions(_autopas, *_particlePropertiesLibrary, _config->deltaT);
      _timers.positionUpdate.stop();
    }
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

    if (_config->deltaT != 0) {
      _timers.velocityUpdate.start();
      TimeDiscretization::calculateVelocities(_autopas, *_particlePropertiesLibrary, _config->deltaT);
      _timers.velocityUpdate.stop();

      // applying Velocity scaling with Thermostat:
      if (_config->useThermostat and (iteration % _config->thermostatInterval) == 0) {
        _timers.thermostat.start();
        Thermostat::apply(_autopas, *_particlePropertiesLibrary, _config->targetTemperature, _config->deltaTemp);
        _timers.thermostat.stop();
      }
    }
  }

  // update temperature for generated config output
  if (_config->useThermostat) {
    _timers.thermostat.start();
    _config->initTemperature = Thermostat::calcTemperature(_autopas, *_particlePropertiesLibrary);
    _timers.thermostat.stop();
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
  cout << fixed << setprecision(_floatStringPrecision);
  cout << endl << "Measurements:" << endl;
  cout << timerToString("Time total      ", durationTotal, digitsTimeTotalMuS);
  cout << timerToString("  Initialization", _timers.init.getTotalTime(), digitsTimeTotalMuS, durationTotal);
  cout << timerToString("  Simulation    ", durationSimulate, digitsTimeTotalMuS, durationTotal);
  cout << timerToString("    Boundaries  ", _timers.boundaries.getTotalTime(), digitsTimeTotalMuS, durationSimulate);
  cout << timerToString("    Position    ", _timers.positionUpdate.getTotalTime(), digitsTimeTotalMuS,
                        durationSimulate);
  cout << timerToString("    Force       ", _timers.forceUpdateTotal.getTotalTime(), digitsTimeTotalMuS,
                        durationSimulate);
  cout << timerToString("      Tuning    ", _timers.forceUpdateTuning.getTotalTime(), digitsTimeTotalMuS,
                        _timers.forceUpdateTotal.getTotalTime());
  cout << timerToString("      NonTuning ", _timers.forceUpdateNonTuning.getTotalTime(), digitsTimeTotalMuS,
                        _timers.forceUpdateTotal.getTotalTime());
  cout << timerToString("    Velocity    ", _timers.velocityUpdate.getTotalTime(), digitsTimeTotalMuS,
                        durationSimulate);
  cout << timerToString("    VTK         ", _timers.vtk.getTotalTime(), digitsTimeTotalMuS, durationSimulate);
  cout << timerToString("    Thermostat  ", _timers.thermostat.getTotalTime(), digitsTimeTotalMuS, durationSimulate);

  auto numIterations = _config->iterations;

  cout << timerToString("One iteration   ", _timers.simulate.getTotalTime() / numIterations, digitsTimeTotalMuS,
                        durationTotal);
  auto mfups = _autopas.getNumberOfParticles(autopas::IteratorBehavior::ownedOnly) * numIterations /
               _timers.forceUpdateTotal.getTotalTime() * 1e-6;
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

std::string Simulation<Particle, ParticleCell>::timerToString(const std::string &name, long timeMS, size_t numberWidth,
                                                              long maxTime) {
  // only print timers that were actually used
  if (timeMS == 0) {
    return "";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(_floatStringPrecision);
  ss << name << " : " << std::setw(numberWidth) << std::right << timeMS << " \u03bcs (" << ((double)timeMS * 1e-6)
     << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeMS / (double)maxTime * 100) << "%";
  }
  ss << std::endl;
  return ss.str();
}
