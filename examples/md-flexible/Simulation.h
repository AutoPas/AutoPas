/**
 * @file Simulation.h
 * @author N. Fottner
 * @date 12/05/19
 */
#pragma once
#include <autopas/utils/MemoryProfiler.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "BoundaryConditions.h"
#include "Generator.h"
#include "PrintableMolecule.h"
#include "Thermostat.h"
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "parsing/MDFlexConfig.h"

template <class Particle, class ParticleCell>
class Simulation {
 public:
  /**
   * Constructor.
   *
   * This starts the timer for total simulation time.
   */
  explicit Simulation() { _timers.startTotal = std::chrono::high_resolution_clock::now(); };

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
    // iterate only over owned Particles, otherwise Simulation explodes
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
  }

  /**Initialized the ParticlePropertiesLibrary for usage in functor
   * during Simualtion::initialize call
   * with the parsed values in the yamlParser
   * */
  void initializeParticlePropertiesLibrary();

  /** @brief This function
   * -initializes the autopas Object with all member speizified in the YamlParser
   * -initializes the simulation domain with the Object Generators
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

  /**Getter for ParticlePropertiesLibrary of Simulation
   * @return unique_prt(ParticlePropertiesLibrary)
   * */
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
    long durationPositionUpdate = 0, durationForceUpdate = 0, durationVelocityUpdate = 0, durationSimulate = 0;
    std::chrono::system_clock::time_point startTotal, stopTotal;
  } _timers;
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
  for (auto &eps : _config->epsilonMap) {
    _particlePropertiesLibrary->addType(eps.first, eps.second, _config->sigmaMap.at(eps.first),
                                        _config->massMap.at(eps.first));
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(const MDFlexConfig &mdFlexConfig) {
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
  _autopas.setCutoff(cutoff);
  _autopas.setVerletSkin(verletSkinRadius);
  _autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  _autopas.setTuningInterval(tuningInterval);
  _autopas.setNumSamples(tuningSamples);
  _autopas.setSelectorStrategy(selectorStrategy);
  _autopas.setAllowedContainers(containerChoice);
  _autopas.setAllowedTraversals(traversalOptions);
  _autopas.setAllowedDataLayouts(dataLayoutOptions);
  _autopas.setAllowedNewton3Options(newton3Options);
  _autopas.setTuningStrategyOption(tuningStrategy);
  _autopas.setAllowedCellSizeFactors(*cellSizeFactors);
  _autopas.setVerletRebuildFrequency(_config->verletRebuildFrequency);
  autopas::Logger::get()->set_level(logLevel);
  _autopas.setBoxMax(_config->boxMax);
  _autopas.setBoxMin(_config->boxMin);
  _autopas.init();
  size_t particleIDCounter = 0;

  // initializing Objects
  for (const auto &grid : cubesGrid) {
    Generator::CubeGrid<Particle, ParticleCell>(_autopas, grid.getTypeId(), particleIDCounter, grid.getBoxMin(),
                                                grid.getParticlesPerDim(), grid.getParticleSpacing(),
                                                grid.getVelocity());
    particleIDCounter += grid.getParticlesTotal();
  }
  for (const auto &cube : cubesGauss) {
    Generator::CubeGauss<Particle, ParticleCell>(_autopas, cube.getTypeId(), particleIDCounter, cube.getBoxMin(),
                                                 cube.getBoxMax(), cube.getParticlesTotal(), cube.getDistributionMean(),
                                                 cube.getDistributionStdDev(), cube.getVelocity());
    particleIDCounter += cube.getParticlesTotal();
  }
  for (const auto &cube : cubesUniform) {
    Generator::CubeRandom<Particle, ParticleCell>(_autopas, cube.getTypeId(), particleIDCounter, cube.getBoxMin(),
                                                  cube.getBoxMax(), cube.getParticlesTotal(), cube.getVelocity());
    particleIDCounter += cube.getParticlesTotal();
  }
  for (const auto &sphere : spheres) {
    Generator::Sphere<Particle, ParticleCell>(_autopas, sphere.getCenter(), sphere.getRadius(),
                                              sphere.getParticleSpacing(), particleIDCounter, sphere.getTypeId(),
                                              sphere.getVelocity());
    particleIDCounter += sphere.getParticlesTotal();
  }

  // initilizing Thermostat
  if (_config->useThermostat) {
    _thermostat = std::make_unique<
        Thermostat<decltype(_autopas), std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>>(
        _config->initTemperature, _config->useCurrentTempForBrownianMotion, _config->targetTemperature,
        _config->deltaTemp, *_particlePropertiesLibrary);
  }

  // initializing velocites of Particles
  if (_config->useThermostat) {
    _thermostat->initialize(_autopas);
  }
}

template <class Particle, class ParticleCell>
template <class FunctorType>
void Simulation<Particle, ParticleCell>::calculateForces() {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
  auto functor = FunctorType(_autopas.getCutoff(), 0.0, *_particlePropertiesLibrary);
  //@todo only iterate over owned particles, right?
  _autopas.iteratePairwise(&functor);
  stopCalc = std::chrono::high_resolution_clock::now();
  auto durationCalcF = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  _timers.durationForceUpdate += durationCalcF;
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate() {
  std::chrono::high_resolution_clock::time_point startSim, stopSim;
  startSim = std::chrono::high_resolution_clock::now();
  // writes initial state of simulation as vtkFile if filename is specified
  if ((not _config->vtkFileName.empty())) {
    this->writeVTKFile(0);
  }
  // main simulation loop
  for (size_t iteration = 0; iteration < _config->iterations; ++iteration) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Iteration " << iteration << std::endl;
    }
    if (_config->periodic) {
      BoundaryConditions<Particle, ParticleCell>::applyPeriodic(_autopas);
    }
    _timers.durationPositionUpdate += _timeDiscretization->CalculateX(_autopas);

    switch (this->_config->functorOption) {
      case MDFlexConfig::FunctorOption::lj12_6: {
        this->calculateForces<autopas::LJFunctor<Particle, ParticleCell, /* mixing */ true>>();
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
    _timers.durationVelocityUpdate += _timeDiscretization->CalculateV(_autopas);
    // applying Velocity scaling with Thermostat:
    if (_config->useThermostat and (iteration % _config->thermostatInterval) == 0) {
      _thermostat->apply(_autopas);
    }
    // only write vtk files periodically and if a filename is given
    if ((not _config->vtkFileName.empty()) and iteration % _config->vtkWriteFrequency == 0) {
      this->writeVTKFile(iteration + 1);
    }
  }

  stopSim = std::chrono::high_resolution_clock::now();
  _timers.durationSimulate = std::chrono::duration_cast<std::chrono::microseconds>(stopSim - startSim).count();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::printStatistics() {
  using namespace std;
  size_t flopsPerKernelCall;

  switch (_config->functorOption) {
    case MDFlexConfig::FunctorOption ::lj12_6: {
      flopsPerKernelCall = autopas::LJFunctor<PrintableMolecule,
                                              autopas::FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_AVX: {
      flopsPerKernelCall =
          autopas::LJFunctorAVX<PrintableMolecule,
                                autopas::FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    default:
      throw std::runtime_error("Invalid Functor choice");
  }

  _timers.stopTotal = std::chrono::high_resolution_clock::now();
  auto durationTotal =
      std::chrono::duration_cast<std::chrono::microseconds>(_timers.stopTotal - _timers.startTotal).count();
  auto durationTotalSec = durationTotal * 1e-6;
  auto durationSimulateSec = _timers.durationSimulate * 1e-6;

  // time statistics
  cout << "Simulation duration without initilization: " << _timers.durationSimulate << " \u03bcs" << endl;
  // Statistics
  cout << fixed << setprecision(3);
  cout << endl << "Measurements:" << endl;
  cout << "Time total   : " << durationTotal << " \u03bcs (" << durationTotalSec << "s)" << endl;
  cout << "Duration of Physics Calculations: " << endl;
  cout << "Force:   " << _timers.durationForceUpdate << " \u03bcs (" << _timers.durationForceUpdate * 1e-6 << "s)"
       << endl;
  cout << "Postion: " << _timers.durationPositionUpdate << " \u03bcs (" << _timers.durationPositionUpdate * 1e-6 << "s)"
       << endl;
  cout << "Velocity " << _timers.durationVelocityUpdate << " \u03bcs (" << _timers.durationVelocityUpdate * 1e-6 << "s)"
       << endl;

  auto numIterations = _config->iterations;

  if (numIterations > 0) {
    cout << "One iteration: " << _timers.durationSimulate / numIterations << " \u03bcs ("
         << durationSimulateSec / numIterations << "s)" << endl;
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
