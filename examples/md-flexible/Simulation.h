//
// Created by nicola on 12.05.19.
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
#include "MDFlexParser.h"
#include "PrintableMolecule.h"
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"

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
  void writeVTKFile(unsigned int iteration) {
    std::string filename = "VtkOutput";
    std::stringstream strstr;
    auto maxNumDigits = std::to_string(_parser->getIterations()).length();
    strstr << filename << "_" << std::setfill('0') << std::setw(maxNumDigits) << iteration << ".vtu";
    std::ofstream vtkFile;
    vtkFile.open(strstr.str());

    vtkFile << "# vtk DataFile Version 2.0" << std::endl;
    vtkFile << "Timestep" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET STRUCTURED_GRID" << std::endl;
    vtkFile << "DIMENSIONS 1 1 1" << std::endl;
    vtkFile << "POINTS " << _autopas.getNumberOfParticles() << " double" << std::endl;

    for (auto iter = _autopas.begin(); iter.isValid(); ++iter) {
      auto pos = iter->getR();
      vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    }

    vtkFile.close();
  }

  /**
   * Constructs a container and fills it with particles.
   *
   * According to the options passed, a DirectSum or 'LinkedCells' container is
   * built. It consists of FullParticleCells and is filled with
   * `PrintableMolecules`. The particles are aligned on a cuboid grid.
   *
   * @param autopas AutoPas object that should be initialized
   * @param particlesPerDim Number of desired particles per dimension.
   * @param particelSpacing Space between two particles along each axis of space.
   */
  void initContainerGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t particlesPerDim,
                         double particelSpacing);

  void initContainerGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, double boxLength, size_t numParticles,
                          double distributionMean, double distributionStdDev);

  void initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas, std::array<double, 3> boxLength,
                            size_t numParticles);

  /**  This function
   * -initializes the autopas Object
   * -sets/initializes the simulation domain with the particles generators
   * @todo -initialized Velocities and Positions (and forces?)
   */
  void initialize(std::shared_ptr<MDFlexParser> parser);

  /**
   * Does the ForceCalculation
   * @param Force Calculation Functor
   * */
  void calculateForces();

  /**
   * This function processes the main simulation loop
   * -calls the time discretization class(calculate fores, etc ...)
   * -do the output each timestep
   * -collects the duration of every Calculation(Position,Force,Velocity)
   */
  void simulate();

  /**Getter for Autopas Oject
   * @return Autopas Object
   */
  autopas::AutoPas<Particle, ParticleCell> *getAutopas() const;
  /**Return current number of Particles in AutoPas Object
   * */
  size_t getNumParticles() { return _autopas.getNumberOfParticles(); }
  /**Prints Statistics(duration of calculation, etc ..) of the Simulation
   * */
  void printStatistics();

 private:
  autopas::AutoPas<Particle, ParticleCell> _autopas;
  std::shared_ptr<MDFlexParser> _parser;
  std::ofstream _logFile;
  std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> _particlePropertiesLibrary;
  std::unique_ptr<TimeDiscretization<decltype(_autopas), double, size_t>> _timeDiscretization;

  struct timers {
    long durationPositionUpdate = 0, durationForceUpdate = 0, durationVelocityUpdate = 0, durationSimulate = 0;
    std::chrono::system_clock::time_point startTotal, stopTotal;
  } _timers;
};

template <class Particle, class ParticleCell>
autopas::AutoPas<Particle, ParticleCell> *Simulation<Particle, ParticleCell>::getAutopas() const {
  return _autopas.get();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(std::shared_ptr<MDFlexParser> parser) {
  _parser = parser;
  double epsilon = _parser->getEpsilon();
  double sigma = _parser->getSigma();
  double mass = _parser->getMass();
  // initialisierung of PCL
  // this implementation doesnt support multiple particle Types, will be coming with PR md-parser
  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibrary>(epsilon, sigma, mass);
  _timeDiscretization =
      std::make_unique<TimeDiscretization<decltype(_autopas)>>(_parser->getDeltaT(), *_particlePropertiesLibrary);

  auto logFileName(_parser->getLogFileName());
  auto particlesTotal(_parser->getParticlesTotal());
  auto particlesPerDim(_parser->getParticlesPerDim());
  auto verletRebuildFrequency(_parser->getVerletRebuildFrequency());
  auto logLevel(_parser->getLogLevel());
  auto &cellSizeFactors(_parser->getCellSizeFactors());
  auto tuningStrategy(_parser->getTuningStrategyOption());
  auto boxLength(_parser->getBoxLength());
  auto containerChoice(_parser->getContainerOptions());
  auto selectorStrategy(_parser->getSelectorStrategy());
  auto cutoff(_parser->getCutoff());
  auto dataLayoutOptions(_parser->getDataLayoutOptions());
  auto distributionMean(_parser->getDistributionMean());
  auto distributionStdDev(_parser->getDistributionStdDev());
  // auto functorChoice(_parser->getFunctorOption());
  auto generatorChoice(_parser->getGeneratorOption());
  auto newton3Options(_parser->getNewton3Options());
  auto particleSpacing(_parser->getParticleSpacing());
  auto traversalOptions(_parser->getTraversalOptions());
  auto tuningInterval(_parser->getTuningInterval());
  auto tuningSamples(_parser->getTuningSamples());
  auto verletSkinRadius(_parser->getVerletSkinRadius());

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

  switch (generatorChoice) {
    case MDFlexParser::GeneratorOption::grid: {
      this->initContainerGrid(_autopas, particlesPerDim,
                              particleSpacing);  // particlesTotal wird in diesem fall in der main geupdated
      break;
    }
    case MDFlexParser::GeneratorOption::uniform: {
      this->initContainerUniform(_autopas, {boxLength, boxLength, boxLength}, particlesTotal);
      break;
    }
    case MDFlexParser::GeneratorOption::gaussian: {
      this->initContainerGauss(_autopas, boxLength, particlesTotal, distributionMean, distributionStdDev);
      break;
    }
    default:
      throw std::runtime_error("Unknown generator choice");
  }
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerGrid(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                           size_t particlesPerDim, double particelSpacing) {
  double ppDxpS = (particlesPerDim)*particelSpacing;
  std::array<double, 3> boxMin({0., 0., 0.});

  std::array<double, 3> boxMax{ppDxpS, ppDxpS, ppDxpS};

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                   {particelSpacing, particelSpacing, particelSpacing},
                                   {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerGauss(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                            double boxLength, size_t numParticles,
                                                            double distributionMean, double distributionStdDev) {
  std::array<double, 3> boxMin({0., 0., 0.});

  std::array<double, 3> boxMax{boxLength, boxLength, boxLength};

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                              std::array<double, 3> boxLength, size_t numParticles) {
  std::array<double, 3> boxMin({0., 0., 0.});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxLength);

  autopas.init();

  PrintableMolecule dummyParticle;
  RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::calculateForces() {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
  auto functor = autopas::LJFunctor<Particle, ParticleCell, /* mixing */ true>(_autopas.getCutoff(), 0.0,
                                                                               *_particlePropertiesLibrary);
  _autopas.iteratePairwise(&functor);
  stopCalc = std::chrono::high_resolution_clock::now();
  auto durationCalcF = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  _timers.durationForceUpdate += durationCalcF;
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate() {
  std::chrono::high_resolution_clock::time_point startSim, stopSim;
  startSim = std::chrono::high_resolution_clock::now();

  // main simulation loop
  for (size_t iteration = 0; iteration < _parser->getIterations(); ++iteration) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Iteration " << iteration << std::endl;
    }

    _timers.durationPositionUpdate += _timeDiscretization.VSCalculateX(_autopas);
    this->calculateForces();

    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << std::endl;
    }
    _timers.durationVelocityUpdate += _timeDiscretization.VSCalculateV(_autopas);

    // only write vtk files periodically and if a filename is given
    if ((not _parser->getVTKFilenName().empty()) and iteration % _parser->getVtkWriteFrequency() == 0) {
      this->writeVTKFile(iteration);
    }
  }

  stopSim = std::chrono::high_resolution_clock::now();
  _timers.durationSimulate = std::chrono::duration_cast<std::chrono::microseconds>(stopSim - startSim).count();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::printStatistics() {
  using namespace std;
  size_t flopsPerKernelCall;

  // FlopsPerKernelCall ließt vom Functor
  switch (_parser->getFunctorOption()) {
    case MDFlexParser::FunctorOption ::lj12_6: {
      flopsPerKernelCall = autopas::LJFunctor<PrintableMolecule,
                                              autopas::FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexParser::FunctorOption ::lj12_6_AVX: {
      flopsPerKernelCall =
          autopas::LJFunctorAVX<PrintableMolecule,
                                autopas::FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    default:
      throw std::runtime_error("Not allowed Functor choice");
  }

  _timers.stopTotal = std::chrono::high_resolution_clock::now();
  auto durationTotal =
      std::chrono::duration_cast<std::chrono::microseconds>(_timers.stopTotal - _timers.startTotal).count();
  auto durationTotalSec = durationTotal * 1e-6;
  auto durationSimulateSec = _timers.durationSimulate * 1e-6;

  // time statistics
  cout << "Simulation duration without initilization: " << _timers.durationSimulate << " \u03bcs" << endl;
  // Statistics
  cout << fixed << setprecision(2);
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
  auto mfups = _autopas.getNumberOfParticles() * numIterations / durationSimulateSec;
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
