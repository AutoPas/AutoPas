//
// Created by nicola on 12.05.19.
//
#pragma once
#include <autopas/utils/MemoryProfiler.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "MDFlexParser.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"

using namespace autopas;
using namespace std;

template <class Particle, class ParticleCell>
class Simulation {
 private:
  shared_ptr<AutoPas<Particle, ParticleCell>> _autopas;
  MDFlexParser *_parser;
  std::ofstream _logFile;

  struct timers {
    long durationPositionUpdate = 0, durationForceUpdate = 0, durationVelocityUpdate = 0, durationSimulate = 0;
    std::chrono::system_clock::time_point startTotal, stopTotal;
  } _timers;

 public:
  explicit Simulation() { _timers.startTotal = std::chrono::high_resolution_clock::now(); };

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
  template <class AutoPasTemplate>
  void writeVTKFile(int iteration, size_t numParticles, AutoPasTemplate &autopas) {
    string filename = "VtkOutput";
    stringstream strstr;
    strstr << filename << "_" << setfill('0') << setw(4) << iteration << ".vtu";
    // string path = "./vtk";
    std::ofstream vtkFile;
    vtkFile.open(strstr.str());

    vtkFile << "# vtk DataFile Version 2.0" << endl;
    vtkFile << "Timestep" << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET STRUCTURED_GRID" << endl;
    vtkFile << "DIMENSIONS 1 1 1" << endl;
    vtkFile << "POINTS " << numParticles << " double" << endl;

    for (auto iter = autopas->begin(); iter.isValid(); ++iter) {
      auto pos = iter->getR();
      vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
    }

    vtkFile.close();
  }

  /**
   * @brief Constructs a container and fills it with particles.
   *
   * According to the options passed, a %DirectSum or %'LinkedCells' container is
   * built. It consists of %`FullParticleCells` and is filled with
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

  void initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas, double boxLength, size_t numParticles);

  /** @brief This function
   * -initializes the autopas Object
   * -sets/initializes the simulation domain with the particles generators
   * @todo -initialized Velocities and Positions (and forces?)
   */
  void initialize(MDFlexParser *parser);

  /**
   * Does the ForceCalculation
   * @param Force Calculation Functor
   * @return Duration of Calculation
   * */
  void CalcF();

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
  AutoPas<Particle, ParticleCell> *getAutopas() const;

  size_t getNumParticles() { return _autopas->getNumberOfParticles(); }

  void printStatistics();
};

template <class Particle, class ParticleCell>
AutoPas<Particle, ParticleCell> *Simulation<Particle, ParticleCell>::getAutopas() const {
  return _autopas.get();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(MDFlexParser *parser) {
  // werte die man später für die initialisierung der Funktoren braucht, temporäre implementierung
  // std::array<double, 3> lowCorner = {0., 0., 0.};
  // std::array<double, 3> highCorner = {5., 5., 5.};
  // double epsilon,sigma  = 1.0;
  _parser = parser;
  auto logFileName(parser->getLogFileName());
  auto particlesTotal(parser->getParticlesTotal());
  auto particlesPerDim(parser->getParticlesPerDim());
  auto verletRebuildFrequency(parser->getVerletRebuildFrequency());
  auto logLevel(parser->getLogLevel());
  auto &cellSizeFactors(parser->getCellSizeFactors());
  auto tuningStrategy(parser->getTuningStrategyOption());
  auto boxLength(parser->getBoxLength());
  auto containerChoice(parser->getContainerOptions());
  auto selectorStrategy(parser->getSelectorStrategy());
  auto cutoff(parser->getCutoff());
  auto dataLayoutOptions(parser->getDataLayoutOptions());
  auto distributionMean(parser->getDistributionMean());
  auto distributionStdDev(parser->getDistributionStdDev());
  // auto functorChoice(parser->getFunctorOption());
  auto generatorChoice(parser->getGeneratorOption());
  auto newton3Options(parser->getNewton3Options());
  auto particleSpacing(parser->getParticleSpacing());
  auto traversalOptions(parser->getTraversalOptions());
  auto tuningInterval(parser->getTuningInterval());
  auto tuningSamples(parser->getTuningSamples());
  auto verletSkinRadius(parser->getVerletSkinRadius());

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
  PrintableMolecule::setEpsilon(parser->getEpsilon());
  PrintableMolecule::setSigma(parser->getSigma());
  PrintableMolecule::setMass(parser->getMass());

  _autopas = make_shared<autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>>>(outputStream);
  _autopas->setCutoff(cutoff);
  _autopas->setVerletSkin(verletSkinRadius);
  _autopas->setVerletRebuildFrequency(verletRebuildFrequency);
  _autopas->setTuningInterval(tuningInterval);
  _autopas->setNumSamples(tuningSamples);
  _autopas->setSelectorStrategy(selectorStrategy);
  _autopas->setAllowedContainers(containerChoice);
  _autopas->setAllowedTraversals(traversalOptions);
  _autopas->setAllowedDataLayouts(dataLayoutOptions);
  _autopas->setAllowedNewton3Options(newton3Options);
  //so that a test in jenkins passes:
  _autopas->setBoxMax({5.,5.,5.});
  //@todo übernommen vom merge: -> prüfen
  _autopas->setTuningStrategyOption(tuningStrategy);
  _autopas->setAllowedCellSizeFactors(cellSizeFactors);
  autopas::Logger::get()->set_level(logLevel);

  switch (generatorChoice) {
    case MDFlexParser::GeneratorOption::grid: {
      this->initContainerGrid(*_autopas, particlesPerDim,
                              particleSpacing);  // particlesTotal wird in diesem fall in der main geupdated
      break;
    }
    case MDFlexParser::GeneratorOption::uniform: {
      this->initContainerUniform(*_autopas, boxLength, particlesTotal);
      break;
    }
    case MDFlexParser::GeneratorOption::gaussian: {
      this->initContainerGauss(*_autopas, boxLength, particlesTotal, distributionMean, distributionStdDev);
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
  std::array<double, 3> boxMax({ppDxpS, ppDxpS, ppDxpS});

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
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                              double boxLength, size_t numParticles) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::CalcF() {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
  //@ TODO: switch for other functors --> mit boolean object?
  //_autopas->iteratePairwise(dynamic_cast<LJFunctor<Particle, ParticleCell>*>(this->_Functor));
  //_autopas->iteratePairwise(this->_Functor);

  auto *functor = new autopas::LJFunctor<Particle, ParticleCell, autopas::FunctorN3Modes::Both, true /* globals */,
                                         true /* relevant for tuning */>(
      _autopas->getCutoff(), 1., 1.0, 0.0);  // cutoff, espi, sigma, shift, box min, box max

  _autopas->iteratePairwise(functor);
  stopCalc = std::chrono::high_resolution_clock::now();
  auto durationCalcF = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  _timers.durationForceUpdate += durationCalcF;
  delete functor;
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate() {
  std::chrono::high_resolution_clock::time_point startSim, stopSim;
  startSim = std::chrono::high_resolution_clock::now();
  double deltaT = _parser->getDeltaT();
  double simTimeNow = 0;
  double simTimeEnd = _parser->getDeltaT() * _parser->getIterations();
  TimeDiscretization<decltype(_autopas)> timeDiscretization(deltaT);

  // main simulation loop
  while (simTimeNow < simTimeEnd) {
    _timers.durationPositionUpdate += timeDiscretization.VSCalculateX(_autopas);

    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      cout << "Iteration " << simTimeNow / deltaT << endl;
      cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
    }
    this->CalcF();
    _timers.durationVelocityUpdate += timeDiscretization.VSCalculateV(_autopas);
    simTimeNow += deltaT;
    this->writeVTKFile(simTimeNow / deltaT, _autopas->getNumberOfParticles(), _autopas);
  }

  stopSim = std::chrono::high_resolution_clock::now();
  _timers.durationSimulate = std::chrono::duration_cast<std::chrono::microseconds>(stopSim - startSim).count();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::printStatistics() {
  size_t flopsPerKernelCall;

  // FlopsPerKernelCall ließt vom Functor
  switch (_parser->getFunctorOption()) {
    case MDFlexParser::FunctorOption ::lj12_6: {
      flopsPerKernelCall =
          LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexParser::FunctorOption ::lj12_6_AVX: {
      flopsPerKernelCall =
          LJFunctorAVX<PrintableMolecule, FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
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

  //@todo conditionally initialize with iterations from simulation time
  auto numIterations = _parser->getIterations();

  if (numIterations > 0) {
    cout << "One iteration: " << _timers.durationSimulate / numIterations << " \u03bcs ("
         << durationSimulateSec / numIterations << "s)" << endl;
  }
  auto mfups = _autopas->getNumberOfParticles() * numIterations / durationSimulateSec;
  cout << "MFUPs/sec    : " << mfups << endl;

  if (_parser->getMeasureFlops()) {
    FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> flopCounterFunctor(
        _autopas->getCutoff());
    _autopas->iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * numIterations;
    // approximation for flops of verlet list generation
    if (_autopas->getContainerType() == autopas::ContainerOption::verletLists)
      flops +=
          flopCounterFunctor.getDistanceCalculations() *
          FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::numFlopsPerDistanceCalculation *
          floor(numIterations / _parser->getVerletRebuildFrequency());

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationSimulateSec << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }
}
