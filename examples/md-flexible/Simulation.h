/**
 * @file Simulation.h
 * @author N. Fottner
 * @date 7/3/19
 */
#pragma once
#include <autopas/utils/MemoryProfiler.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "VTKWriter/VTKWriter.h"
#include "BoundaryConditions.h"
#include "Generator.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "TimeDiscretization.h"
#include "YamlParser.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"

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
    if (not _parser->getLogFileName().empty()) {
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
  /**Writes a VTK file for the current state of the AutoPas object
   * code taken from the MolSim course
   * @param iteration
   * */
  void plotVTK(int iteration) {
      writer.initializeOutput(_autopas.getNumberOfParticles());
      for(auto iter = _autopas.begin();iter.isValid();++iter){
          writer.plotParticle(*iter);
      }
      writer.writeFile(_parser->getVTKFileName(), iteration);
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
  void initialize(const std::shared_ptr<YamlParser> &parser);

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
   * Return the current number of Particles in the AutoPas Object.
   * @return
   */
  size_t getNumParticles() { return _autopas.getNumberOfParticles(); }

  /**
   * Prints statistics like duration of calculation etc of the Simulation.
   */
  void printStatistics();
  /**Setter for parser Object
   * @param parser
   * */
  void setParser(const YamlParser &parser);
  /**Getter for ParticlePropertiesLibrary of Simulation
   * @return unique_prt(ParticlePropertiesLibrary)
   * */
  const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &getPpl() const;

 private:
  outputWriter::VTKWriter writer;
  autopas::AutoPas<Particle, ParticleCell> _autopas;
  std::shared_ptr<YamlParser> _parser;
  std::ofstream _logFile;
  std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> _particlePropertiesLibrary;
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
  _particlePropertiesLibrary = std::make_unique<std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>();
  std::map<unsigned long, double> epsilonMap = _parser->getEpsilonMap();
  std::map<unsigned long, double> sigmaMap = _parser->getSigmaMap();
  std::map<unsigned long, double> massMap = _parser->getMassMap();
  if (epsilonMap.empty()) {
    // initializing PPL with default values epsilon=sigma=mass=1.0
    double epsilon = 1.0;
    double sigma = 1.0;
    double mass = 1.0;
    _particlePropertiesLibrary->addType(0, epsilon, sigma, mass);
  } else {
    // all 3 maps are well initialized in parser(parser catches error and throws exception if something is going wrong)
    for (auto eps : epsilonMap) {
      _particlePropertiesLibrary->addType(eps.first, eps.second, sigmaMap.at(eps.first), massMap.at(eps.first));
    }
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(const std::shared_ptr<YamlParser> &parser) {
  _parser = parser;
  initializeParticlePropertiesLibrary();
  _timeDiscretization = std::make_unique<
      TimeDiscretization<decltype(_autopas), std::remove_reference_t<decltype(*_particlePropertiesLibrary)>>>(
      _parser->getDeltaT(), *_particlePropertiesLibrary);
  auto logFileName(_parser->getLogFileName());
  auto verletRebuildFrequency(_parser->getVerletRebuildFrequency());
  auto logLevel(_parser->getLogLevel());
  auto &cellSizeFactors(_parser->getCellSizeFactors());
  auto tuningStrategy(_parser->getTuningStrategyOption());
  auto containerChoice(_parser->getContainerOptions());
  auto selectorStrategy(_parser->getSelectorStrategy());
  auto cutoff(_parser->getCutoff());
  auto dataLayoutOptions(_parser->getDataLayoutOptions());
  auto newton3Options(_parser->getNewton3Options());
  auto traversalOptions(_parser->getTraversalOptions());
  auto tuningInterval(_parser->getTuningInterval());
  auto tuningSamples(_parser->getTuningSamples());
  auto verletSkinRadius(_parser->getVerletSkinRadius());
  auto CubeGrid(_parser->getCubeGrid());
  auto CubeGauss(_parser->getCubeGauss());
  auto CubeUniform(_parser->getCubeUniform());
  auto Sphere(_parser->getSphere());
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
  _autopas.setAllowedCellSizeFactors(cellSizeFactors);
  autopas::Logger::get()->set_level(logLevel);
  _autopas.setBoxMax(_parser->getBoxMax());
  _autopas.setBoxMin(_parser->getBoxMin());
  _autopas.init();
  size_t idcounter = 0;
  for (auto C : CubeGrid) {
    Generator::CubeGrid<Particle, ParticleCell>(_autopas, C.getTypeId(), idcounter, C.getBoxMin(),
                                                C.getParticlesPerDim(), C.getParticleSpacing(), C.getVelocity());
    idcounter = +C.getParticlesTotal();
  }
  for (auto C : CubeGauss) {
    Generator::CubeGauss<Particle, ParticleCell>(_autopas, C.getTypeId(), idcounter, C.getBoxMin(), C.getBoxMax(),
                                                 C.getParticlesTotal(), C.getDistributionMean(),
                                                 C.getDistributionStdDev(), C.getVelocity());
    idcounter = +C.getParticlesTotal();
  }
  for (auto C : CubeUniform) {
    Generator::CubeRandom<Particle, ParticleCell>(_autopas, C.getTypeId(), idcounter, C.getBoxMin(), C.getBoxMax(),
                                                  C.getParticlesTotal(), C.getVelocity());
    idcounter = +C.getParticlesTotal();
  }
  for (auto S : Sphere) {
    Generator::Sphere<Particle, ParticleCell>(_autopas, S.getCenter(), S.getRadius(), S.getParticleSpacing(), idcounter,
                                              S.getTypeId(), S.getVelocity());
    idcounter = +S.getParticlesTotal();
  }
  for(auto iter=_autopas.begin();iter.isValid();++iter){
      std::cout <<iter->toString() << std::endl;
  }
}

template <class Particle, class ParticleCell>
template <class FunctorType>
void Simulation<Particle, ParticleCell>::calculateForces() {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
  auto functor = FunctorType(_autopas.getCutoff(), 0.0, *_particlePropertiesLibrary);
  _autopas.iteratePairwise(&functor);
  stopCalc = std::chrono::high_resolution_clock::now();
  auto durationCalcF = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  _timers.durationForceUpdate += durationCalcF;
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate() {
  std::chrono::high_resolution_clock::time_point startSim, stopSim;
  startSim = std::chrono::high_resolution_clock::now();
  //writes initial state of simulation as vtkFile if filename is specified
    if ((not _parser->getVTKFileName().empty())) {
        this->plotVTK(0);
    }
    // main simulation loop
  for (size_t iteration = 0; iteration < _parser->getIterations(); ++iteration) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Iteration " << iteration << std::endl;
    }
      if(_parser->isPeriodic()) {
          BoundaryConditions<Particle,ParticleCell>::applyPeriodic(_autopas);
      }
    _timers.durationPositionUpdate += _timeDiscretization->CalculateX(_autopas);

    switch (this->_parser->getFunctorOption()) {
      case YamlParser::FunctorOption::lj12_6: {
        this->calculateForces<autopas::LJFunctor<Particle, ParticleCell, /* mixing */ true>>();
        break;
      }
      case YamlParser::FunctorOption::lj12_6_AVX: {
        this->calculateForces<autopas::LJFunctorAVX<Particle, ParticleCell, /* mixing */ true>>();
        break;
      }
    }
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << std::endl;
    }
    _timers.durationVelocityUpdate += _timeDiscretization->CalculateV(_autopas);

    // only write vtk files periodically and if a filename is given
    if ((not _parser->getVTKFileName().empty()) and iteration % _parser->getVtkWriteFrequency() == 0) {
      this->plotVTK(iteration+1);
    }
  }

  stopSim = std::chrono::high_resolution_clock::now();
  _timers.durationSimulate = std::chrono::duration_cast<std::chrono::microseconds>(stopSim - startSim).count();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::printStatistics() {
  using namespace std;
  size_t flopsPerKernelCall;

  switch (_parser->getFunctorOption()) {
    case YamlParser::FunctorOption ::lj12_6: {
      flopsPerKernelCall = autopas::LJFunctor<PrintableMolecule,
                                              autopas::FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    case YamlParser::FunctorOption ::lj12_6_AVX: {
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

  auto numIterations = _parser->getIterations();

  if (numIterations > 0) {
    cout << "One iteration: " << _timers.durationSimulate / numIterations << " \u03bcs ("
         << durationSimulateSec / numIterations << "s)" << endl;
  }
  auto mfups = _autopas.getNumberOfParticles() * numIterations / durationSimulateSec * 1e-6;
  cout << "MFUPs/sec    : " << mfups << endl;

  if (_parser->getMeasureFlops()) {
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
          floor(numIterations / _parser->getVerletRebuildFrequency());

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationSimulateSec << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::setParser(const YamlParser &parser) {
  _parser = std::make_shared<YamlParser>(parser);
}

template <class Particle, class ParticleCell>
const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &Simulation<Particle, ParticleCell>::getPpl() const {
  return _particlePropertiesLibrary;
}
