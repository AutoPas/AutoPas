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
#include "Generator.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "TimeDiscretization.h"
#include "YamlParser.h"
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"
using namespace autopas;
using namespace std;

template <class Particle, class ParticleCell>
class Simulation {
 private:
  AutoPas<Particle, ParticleCell> _autopas;
  YamlParser _parser;
  std::ofstream _logFile;
  shared_ptr<ParticleClassLibrary> _PCL;

  struct timers {
    long durationPositionUpdate = 0, durationForceUpdate = 0, durationVelocityUpdate = 0, durationSimulate = 0;
    std::chrono::system_clock::time_point startTotal, stopTotal;
  } _timers;

 public:
  explicit Simulation() { _timers.startTotal = std::chrono::high_resolution_clock::now(); };

  ~Simulation() {
    if (not _parser.getLogFileName().empty()) {
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
  void writeVTKFile(int iteration, size_t numParticles, autopas::AutoPas<Particle,ParticleCell>  &autopas) {
    string filename = "VtkOutput";
    stringstream strstr;
    strstr << filename << "_" << setfill('0') << setw(4) << iteration << ".vtu";
    // string path = "./vtk";
    std::ofstream vtkFile;
    vtkFile.open(strstr.str());

    vtkFile << "# vtk DataFile Version 2.0" << endl;
    vtkFile << "Timestep" << endl;
    vtkFile << "ASCII" << endl;
    //@todo adapt to current state
    vtkFile << "DATASET STRUCTURED_GRID" << endl;
    vtkFile << "DIMENSIONS 1 1 1" << endl;
    vtkFile << "POINTS " << numParticles << " double" << endl;

    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
      auto pos = iter->getR();
      vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
    }

    vtkFile.close();
  }

  /** @brief This function
   * -initializes the autopas Object
   * -sets/initializes the simulation domain with the particles generators
   * @todo -initialized Velocities and Positions (and forces?)
   */
  void initialize(YamlParser &parser);

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

  size_t getNumParticles() { return _autopas.getNumberOfParticles(); }

  void printStatistics();
};

template <class Particle, class ParticleCell>
AutoPas<Particle, ParticleCell> *Simulation<Particle, ParticleCell>::getAutopas() const {
  return _autopas.get();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(YamlParser &parser) {
  //@todo schöner machen:
  _parser = parser;
  //@todo lösen:
  // temporäre implemetierung mit nur einer particle Class
  double epsilon = _parser.getEpsilon();
  double sigma = _parser.getSigma();
  double mass = _parser.getMass();
  //@todo PCL richtig initialisieren
  _PCL = make_shared<ParticleClassLibrary>(epsilon, sigma, mass, _parser.particlesTotal());
  // initialisierung of
  auto logFileName(_parser.getLogFileName());
  auto particlesTotal(_parser.particlesTotal()) ;
  auto verletRebuildFrequency(_parser.getVerletRebuildFrequency());
  auto logLevel(_parser.getLogLevel());
  auto &cellSizeFactors(_parser.getCellSizeFactors());
  auto tuningStrategy(_parser.getTuningStrategyOption());
  auto containerChoice(_parser.getContainerOptions());
  auto selectorStrategy(_parser.getSelectorStrategy());
  auto cutoff(_parser.getCutoff());
  auto dataLayoutOptions(_parser.getDataLayoutOptions());
  auto newton3Options(_parser.getNewton3Options());
  auto traversalOptions(_parser.getTraversalOptions());
  auto tuningInterval(_parser.getTuningInterval());
  auto tuningSamples(_parser.getTuningSamples());
  auto verletSkinRadius(_parser.getVerletSkinRadius());
  auto CubeGrid(_parser.getCubeGrid());
  auto CubeGauss(_parser.getCubeGauss());
  auto CubeUniform(_parser.getCubeUniform());
  auto Sphere(_parser.getSphere());

  // select either std::out or a logfile for autopas log output.
  // This does not affect md-flex output.
  std::streambuf *streamBuf;
  if (logFileName.empty()) {
    streamBuf = std::cout.rdbuf();
  } else {
    _logFile.open(logFileName);
    streamBuf = _logFile.rdbuf();
  }
  //@todo autopas mit logfilename richtig initialisieren
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
  _autopas.setBoxMax(_parser.getBoxMax());
    _autopas.setBoxMin(_parser.getBoxMin());
    _autopas.init();



    for (auto C : CubeGrid) {
    Generator::CubeGrid<Particle,ParticleCell>(_autopas,C.getBoxMin(), C.getParticlesPerDim(), C.getParticleSpacing(), C.getVelocity());
  }
  for (auto C : CubeGauss) {
    Generator::CubeGauss<Particle,ParticleCell>(_autopas,C.getBoxMin(),C.getBoxMax(), C.getNumParticles(), C.getDistributionMean(),
                         C.getDistributionStdDev(), C.getVelocity());
  }
  for (auto C : CubeUniform) {
    Generator::CubeRandom<Particle,ParticleCell>(_autopas,C.getBoxMin(),C.getBoxMax(), C.getNumParticles(), C.getVelocity());
  }
  for (auto S : Sphere) {
    Generator::Sphere<Particle,ParticleCell>(_autopas, S.getCenter(), S.getRadius(), S.getParticleSpacing(), S.getId(), S.getVelocity());
  }

}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::CalcF() {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
  //@ TODO: switch for other functors --> mit boolean object?
  auto functor = LJFunctor<Particle, ParticleCell>(_autopas.getCutoff(), *_PCL, 0.0);
  _autopas.iteratePairwise(&functor);
  stopCalc = std::chrono::high_resolution_clock::now();
  auto durationCalcF = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  _timers.durationForceUpdate += durationCalcF;
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::simulate() {
  std::chrono::high_resolution_clock::time_point startSim, stopSim;
  startSim = std::chrono::high_resolution_clock::now();
  double deltaT = _parser.getDeltaT();
  double simTimeNow = 0;
  double simTimeEnd = _parser.getDeltaT() * _parser.getIterations();
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
    this->writeVTKFile(simTimeNow / deltaT, _autopas.getNumberOfParticles(), _autopas);
  }

  stopSim = std::chrono::high_resolution_clock::now();
  _timers.durationSimulate = std::chrono::duration_cast<std::chrono::microseconds>(stopSim - startSim).count();
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::printStatistics() {
  size_t flopsPerKernelCall;

  // FlopsPerKernelCall ließt vom Functor
  switch (_parser.getFunctorOption()) {
    case YamlParser::FunctorOption ::lj12_6: {
      flopsPerKernelCall =
          LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    case YamlParser::FunctorOption ::lj12_6_AVX: {
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

  auto numIterations = _parser.getIterations();

  if (numIterations > 0) {
    cout << "One iteration: " << _timers.durationSimulate / numIterations << " \u03bcs ("
         << durationSimulateSec / numIterations << "s)" << endl;
  }
  auto mfups = _autopas.getNumberOfParticles() * numIterations / durationSimulateSec;
  cout << "MFUPs/sec    : " << mfups << endl;

  if (_parser.getMeasureFlops()) {
    FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> flopCounterFunctor(_autopas.getCutoff());
    _autopas.iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * numIterations;
    // approximation for flops of verlet list generation
    if (_autopas.getContainerType() == autopas::ContainerOption::verletLists)
      flops +=
          flopCounterFunctor.getDistanceCalculations() *
          FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::numFlopsPerDistanceCalculation *
          floor(numIterations / _parser.getVerletRebuildFrequency());

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationSimulateSec << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }
}
