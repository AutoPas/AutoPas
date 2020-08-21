/**
 * @file Simulation.h
 * @author N. Fottner
 * @date 12/05/19
 */
#pragma once

#include <fstream>
#include <iostream>

#include "BoundaryConditions.h"
#include "Checkpoint.h"
#include "Generator.h"
#include "PrintableMolecule.h"
#include "Thermostat.h"
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "autopas/utils/MemoryProfiler.h"
#include "autopasTools/generators/GaussianGenerator.h"
#include "autopasTools/generators/GridGenerator.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "parsing/MDFlexConfig.h"

/**
 * The main simulation class.
 * @tparam Particle
 * @tparam ParticleCell
 */
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
   * @param iteration
   * @param autopas
   */
  void writeVTKFile(unsigned int iteration, autopas::AutoPas<Particle, ParticleCell> &autopas);

  /**
   * Initializes the ParticlePropertiesLibrary with properties from _config.
   */
  void initializeParticlePropertiesLibrary();

  /**
   * Initializes the AutoPas Object with the given config and initializes the simulation domain with the Object
   * Generators.
   * @param mdFlexConfig
   * @param autopas
   */
  void initialize(const MDFlexConfig &mdFlexConfig, autopas::AutoPas<Particle, ParticleCell> &autopas);

  /**
   * Calculates the pairwise forces in the system and measures the runtime.
   * @tparam Force Calculation Functor
   * @param autopas
   */
  template <class FunctorType>
  void calculateForces(autopas::AutoPas<Particle, ParticleCell> &autopas);

  /**
   * Calculate influences from global, non pairwise forces, e.g. gravity.
   * @param autopas
   */
  void globalForces(autopas::AutoPas<Particle, ParticleCell> &autopas);

  /**
   * This function processes the main simulation loop
   * -calls the time discretization class(calculate fores, etc ...)
   * -do the output each timestep
   * -collects the duration of every Calculation(Position,Force,Velocity)
   * @param autopas
   */
  void simulate(autopas::AutoPas<Particle, ParticleCell> &autopas);

  /**
   * Indicates if enough iterations were completed yet.
   * Uses class member variables.
   * @return
   */
  [[nodiscard]] bool needsMoreIterations() const;

  /**
   * Prints statistics like duration of calculation etc of the Simulation.
   * @param autopas
   */
  void printStatistics(autopas::AutoPas<Particle, ParticleCell> &autopas);

  /**
   * Getter for ParticlePropertiesLibrary of Simulation.
   * @return unique_prt(ParticlePropertiesLibrary)
   */
  [[nodiscard]] const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &getPpl() const;

  /**
   * Calculate the homogeneity of the scenario by using the standard deviation.
   * @param autopas
   * @return double
   */
  double calculateHomogeneity(autopas::AutoPas<Particle, ParticleCell> &autopas);

 private:
  using ParticlePropertiesLibraryType = ParticlePropertiesLibrary<double, size_t>;
  constexpr static bool _shifting = true;
  constexpr static bool _mixing = true;

  std::shared_ptr<MDFlexConfig> _config;
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;

  struct timers {
    autopas::utils::Timer positionUpdate, forceUpdateTotal, forceUpdatePairwise, forceUpdateGlobal, forceUpdateTuning,
        forceUpdateNonTuning, velocityUpdate, simulate, vtk, init, total, thermostat, boundaries;
  } _timers;

  /**
   * Number of completed iterations. Aka. number of current iteration.
   */
  size_t iteration = 0;
  /**
   * Counts completed iterations that were used for tuning
   */
  size_t numTuningIterations = 0;
  /**
   * Counts completed tuning phases.
   */
  size_t numTuningPhasesCompleted = 0;
  /**
   * Indicator if the previous iteration was used for tuning.
   */
  bool previousIterationWasTuningIteration = false;

  /**
   * Precision of floating point numbers printed.
   */
  constexpr static auto _floatStringPrecision = 3;

  /**
   * Homogeneity of the scenario, calculated by the standard deviation of the density.
   */
  double homogeneity = 0;

  /**
   * Convert a time and a name to a properly formatted string.
   * @param name incl. offset.
   * @param timeNS in nanoseconds.
   * @param numberWidth Width to which the time should be offset.
   * @param maxTime if passed the percentage of timeNS of maxTime is appended.
   * @return formatted std::string
   */
  std::string timerToString(const std::string &name, long timeNS, size_t numberWidth = 0, long maxTime = 0);
};

template <typename Particle, typename ParticleCell>
void Simulation<Particle, ParticleCell>::initializeParticlePropertiesLibrary() {
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
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(const MDFlexConfig &mdFlexConfig,
                                                    autopas::AutoPas<Particle, ParticleCell> &autopas) {
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
  autopas.setEvidenceFirstPrediction(_config->evidenceFirstPrediction.value);
  autopas.setExtrapolationMethodOption(_config->extrapolationMethodOption.value);
  autopas.setNumSamples(_config->tuningSamples.value);
  autopas.setMaxEvidence(_config->tuningMaxEvidence.value);
  autopas.setSelectorStrategy(_config->selectorStrategy.value);
  autopas.setTuningInterval(_config->tuningInterval.value);
  autopas.setTuningStrategyOption(_config->tuningStrategyOption.value);
  autopas.setVerletClusterSize(_config->verletClusterSize.value);
  autopas.setVerletRebuildFrequency(_config->verletRebuildFrequency.value);
  autopas.setVerletSkin(_config->verletSkinRadius.value);
  autopas.setAcquisitionFunction(_config->acquisitionFunctionOption.value);
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
  for (const auto &grid : _config->cubeGridObjects) {
    Generator::cubeGrid<Particle, ParticleCell>(autopas, grid);
  }
  for (const auto &cube : _config->cubeGaussObjects) {
    Generator::cubeGauss<Particle, ParticleCell>(autopas, cube);
  }
  for (const auto &cube : _config->cubeUniformObjects) {
    Generator::cubeRandom<Particle, ParticleCell>(autopas, cube);
  }
  for (const auto &sphere : _config->sphereObjects) {
    Generator::sphere<Particle, ParticleCell>(autopas, sphere);
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

template <class Particle, class ParticleCell>
template <class FunctorType>
void Simulation<Particle, ParticleCell>::calculateForces(autopas::AutoPas<Particle, ParticleCell> &autopas) {
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
    ++numTuningIterations;
  } else {
    _timers.forceUpdateNonTuning.addTime(timeIteration);
    // if the previous iteration was a tuning iteration an the current one is not
    // we have reached the end of a tuning phase
    if (previousIterationWasTuningIteration) {
      ++numTuningPhasesCompleted;
    }
  }
  previousIterationWasTuningIteration = tuningIteration;

  // global forces

  _timers.forceUpdateGlobal.start();
  globalForces(autopas);
  _timers.forceUpdateGlobal.stop();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::globalForces(autopas::AutoPas<Particle, ParticleCell> &autopas) {
  // skip application of zero force
  if (_config->globalForceIsZero()) {
    return;
  }

#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(autopas)
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    iter->setV(autopas::utils::ArrayMath::add(iter->getV(), _config->globalForce.value));
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate(autopas::AutoPas<Particle, ParticleCell> &autopas) {
  this->homogeneity = Simulation::calculateHomogeneity(autopas);
  _timers.simulate.start();

  // main simulation loop
  for (; needsMoreIterations(); ++iteration) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Iteration " << iteration << std::endl;
    }

    if (_config->deltaT.value != 0) {
      // only write vtk files periodically and if a filename is given.
      if ((not _config->vtkFileName.value.empty()) and iteration % _config->vtkWriteFrequency.value == 0) {
        this->writeVTKFile(iteration, autopas);
      }

      // calculate new positions
      _timers.positionUpdate.start();
      TimeDiscretization::calculatePositions(autopas, *_particlePropertiesLibrary, _config->deltaT.value);
      _timers.positionUpdate.stop();

      // apply boundary conditions AFTER the position update!
      if (_config->periodic.value) {
        _timers.boundaries.start();
        BoundaryConditions<ParticleCell>::applyPeriodic(autopas, false);
        _timers.boundaries.stop();
      } else {
        throw std::runtime_error(
            "Simulation::simulate(): at least one boundary condition has to be set. Please enable the periodic "
            "boundary conditions!");
      }
    }
    switch (this->_config->functorOption.value) {
      case MDFlexConfig::FunctorOption::lj12_6: {
        this->calculateForces<autopas::LJFunctor<Particle, ParticleCell, _shifting, _mixing>>(autopas);
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_Globals: {
        this->calculateForces<autopas::LJFunctor<Particle, ParticleCell, _shifting, _mixing,
                                                 autopas::FunctorN3Modes::Both, /* globals */ true>>(autopas);
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_AVX: {
        this->calculateForces<autopas::LJFunctorAVX<Particle, ParticleCell, _shifting, _mixing>>(autopas);
        break;
      }
    }
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << std::endl;
    }

    if (_config->deltaT.value != 0) {
      _timers.velocityUpdate.start();
      TimeDiscretization::calculateVelocities(autopas, *_particlePropertiesLibrary, _config->deltaT.value);
      _timers.velocityUpdate.stop();

      // applying Velocity scaling with Thermostat:
      if (_config->useThermostat.value and (iteration % _config->thermostatInterval.value) == 0) {
        _timers.thermostat.start();
        Thermostat::apply(autopas, *_particlePropertiesLibrary, _config->targetTemperature.value,
                          _config->deltaTemp.value);
        _timers.thermostat.stop();
      }
    }
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
    BoundaryConditions<ParticleCell>::applyPeriodic(autopas, true);
    _timers.boundaries.stop();
    this->writeVTKFile(_config->iterations.value, autopas);
  }

  _timers.simulate.stop();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::printStatistics(autopas::AutoPas<Particle, ParticleCell> &autopas) {
  using namespace std;
  size_t flopsPerKernelCall = 0;

  switch (_config->functorOption.value) {
    case MDFlexConfig::FunctorOption ::lj12_6: {
      flopsPerKernelCall = autopas::LJFunctor<Particle, ParticleCell, _shifting, _mixing>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_Globals: {
      flopsPerKernelCall = autopas::LJFunctor<Particle, ParticleCell, _shifting, _mixing, autopas::FunctorN3Modes::Both,
                                              /* globals */ true>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_AVX: {
      flopsPerKernelCall =
          autopas::LJFunctorAVX<Particle, ParticleCell, _shifting, _mixing>::getNumFlopsPerKernelCall();
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

  cout << timerToString("One iteration   ", _timers.simulate.getTotalTime() / iteration, digitsTimeTotalNS,
                        durationTotal);
  auto mfups = autopas.getNumberOfParticles(autopas::IteratorBehavior::ownedOnly) * iteration * 1e-6 /
               (_timers.forceUpdateTotal.getTotalTime() * 1e-9);  // 1e-9 for ns to s, 1e-6 for M in MFUP
  cout << "Tuning iterations: " << numTuningIterations << " / " << iteration << " = "
       << ((double)numTuningIterations / iteration * 100) << "%" << endl;
  cout << "MFUPs/sec    : " << mfups << endl;
  cout << "Standard Deviation of Homogeneity    : " << homogeneity << endl;

  if (_config->dontMeasureFlops.value) {
    autopas::FlopCounterFunctor<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> flopCounterFunctor(
        autopas.getCutoff());
    autopas.iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * iteration;
    // approximation for flops of verlet list generation
    if (autopas.getContainerType() == autopas::ContainerOption::verletLists)
      flops +=
          flopCounterFunctor.getDistanceCalculations() *
          autopas::FlopCounterFunctor<PrintableMolecule,
                                      autopas::FullParticleCell<PrintableMolecule>>::numFlopsPerDistanceCalculation *
          floor(iteration / _config->verletRebuildFrequency.value);

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

std::string Simulation<Particle, ParticleCell>::timerToString(const std::string &name, long timeNS, size_t numberWidth,
                                                              long maxTime) {
  // only print timers that were actually used
  if (timeNS == 0) {
    return "";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(_floatStringPrecision);
  ss << name << " : " << std::setw(numberWidth) << std::right << timeNS << " ns (" << ((double)timeNS * 1e-9) << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeNS / (double)maxTime * 100) << "%";
  }
  ss << std::endl;
  return ss.str();
}
template <class Particle, class ParticleCell>
bool Simulation<Particle, ParticleCell>::needsMoreIterations() const {
  return iteration < _config->iterations.value or numTuningPhasesCompleted < _config->tuningPhases.value;
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::writeVTKFile(unsigned int iteration,
                                                      autopas::AutoPas<Particle, ParticleCell> &autopas) {
  _timers.vtk.start();

  std::string fileBaseName = _config->vtkFileName.value;
  // only count number of owned particles here
  const auto numParticles = autopas.getNumberOfParticles(autopas::IteratorBehavior::ownedOnly);
  std::ostringstream strstr;
  auto maxNumDigits = std::to_string(_config->iterations.value).length();
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
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  }
  vtkFile << std::endl;

  vtkFile << "POINT_DATA " << numParticles << std::endl;
  // print velocities
  vtkFile << "VECTORS velocities double" << std::endl;
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto v = iter->getV();
    vtkFile << v[0] << " " << v[1] << " " << v[2] << std::endl;
  }
  vtkFile << std::endl;

  // print Forces
  vtkFile << "VECTORS forces double" << std::endl;
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto f = iter->getF();
    vtkFile << f[0] << " " << f[1] << " " << f[2] << std::endl;
  }
  vtkFile << std::endl;

  // print TypeIDs
  vtkFile << "SCALARS typeIds int" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    vtkFile << iter->getTypeId() << std::endl;
  }
  vtkFile << std::endl;

  vtkFile.close();

  _timers.vtk.stop();
}

template <class Particle, class ParticleCell>
double Simulation<Particle, ParticleCell>::calculateHomogeneity(autopas::AutoPas<Particle, ParticleCell> &autopas) {
  int numberOfParticles = autopas.getNumberOfParticles();
  // approximately the resolution we want to get.
  int numberOfCells = ceil(numberOfParticles / 10.);

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
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    std::array<double, 3> particleLocation = iter->getR();
    std::array<size_t, 3> index = {};
    for (int i = 0; i < particleLocation.size(); i++) {
      index[i] = particleLocation[i] / cellLength;
    }
    const size_t cellIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(index, cellsPerDimension);
    particlesPerCell[cellIndex] += 1;
    // calculate the size of the current cell
    allVolumes[cellIndex] = 0;
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
    densityPerCell[i] = (particlesPerCell[i] == 0)
                            ? 0
                            : (particlesPerCell[i] / allVolumes[i]);  // make sure there is no division of zero
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
