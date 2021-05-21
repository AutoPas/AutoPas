/**
 * @file MDFlexSingleRank.h
 * @author J. Körner
 * @date 07.04.2021
 */
#include "MDFlexSingleRank.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "autopas/utils/MemoryProfiler.h"
#include "src/BoundaryConditions.h"
#include "src/Thermostat.h"

MDFlexSingleRank::MDFlexSingleRank(int dimensionCount, int argc, char **argv){
	MDFlexSimulation::initialize(dimensionCount, argc, argv);
}

void MDFlexSingleRank::run() {
  std::cout << std::endl << "Using " << autopas::autopas_get_max_threads() << " Threads" << std::endl;
  std::cout << "Starting simulation... " << std::endl;
	
  _timers.simulate.start();

  auto [maxIterationsEstimate, maxIterationsIsPrecise] = estimateNumberOfIterations();

  // main simulation loop
  for (; needsMoreIterations(); ++_iteration) {
    if (not _configuration->dontShowProgressBar.value) {
      printProgress(_iteration, maxIterationsEstimate, maxIterationsIsPrecise);
    }

    // only do time step related stuff when there actually is time-stepping
    if (_configuration->deltaT.value != 0) {
      // only write vtk files periodically and if a filename is given.
      if (not _configuration->vtkFileName.value.empty() and _iteration % _configuration->vtkWriteFrequency.value == 0) {
        this->writeVTKFile();
      }

      // calculate new positions
      _timers.positionUpdate.start();
      updatePositions();
      _timers.positionUpdate.stop();

      // apply boundary conditions AFTER the position update!
      if (_configuration->periodic.value) {
        _timers.boundaries.start();
        BoundaryConditions::applyPeriodic(*_autoPasContainer, false);
        _timers.boundaries.stop();
      } else {
        throw std::runtime_error(
            "Simulation::simulate(): at least one boundary condition has to be set. Please enable the periodic "
            "boundary conditions!");
      }
    }

		updateForces();

    // only show memory usage in when the logger is set to debug
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << std::endl;
    }

		updateVelocities();
  }

  // final update for a full progress bar
  if (not _configuration->dontShowProgressBar.value) {
    // The last update is precise, so we know the number of iterations.
    printProgress(_iteration, _iteration, true);
    // The progress bar does not end the line. Since this is the last progress bar, end the line here.
    std::cout << std::endl;
  }

  // update temperature for generated config output
  if (_configuration->useThermostat.value) {
    _timers.thermostat.start();
    _configuration->initTemperature.value = Thermostat::calcTemperature(*_autoPasContainer,
			*(_configuration->getParticlePropertiesLibrary()));
    _timers.thermostat.stop();
  }

  // writes final state of the simulation
  if ((not _configuration->vtkFileName.value.empty())) {
    _timers.boundaries.start();
    BoundaryConditions::applyPeriodic(*_autoPasContainer, true);
    _timers.boundaries.stop();
    this->writeVTKFile();
  }

  _timers.simulate.stop();
  std::cout << "Simulation done!" << std::endl << std::endl;

  // Statistics about the simulation
  printStatistics();
}

void MDFlexSingleRank::initializeDomainDecomposition(int &dimensionCount){
	std::vector<double> boxMin(_configuration->boxMin.value.begin(), _configuration->boxMin.value.end());
	std::vector<double> boxMax(_configuration->boxMax.value.begin(), _configuration->boxMax.value.end());
	
	_domainDecomposition = std::make_shared<SingleDomain>(_argc, _argv, dimensionCount, boxMin, boxMax);

	std::vector<double> localBoxMin = _domainDecomposition->getLocalBoxMin();
	std::vector<double> localBoxMax = _domainDecomposition->getLocalBoxMax();
	
	for (int i = 0; i < localBoxMin.size(); ++i){
		_configuration->boxMin.value[i] = localBoxMin[i];
		_configuration->boxMax.value[i] = localBoxMax[i];
	}
}

void MDFlexSingleRank::printStatistics() {
  using namespace std;
  size_t flopsPerKernelCall;

  switch (_configuration->functorOption.value) {
    case MDFlexConfig::FunctorOption ::lj12_6: {
      flopsPerKernelCall = autopas::LJFunctor<ParticleType, true, true>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_Globals: {
      flopsPerKernelCall = autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                              /* globals */ true>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexConfig::FunctorOption ::lj12_6_AVX: {
      flopsPerKernelCall = autopas::LJFunctorAVX<ParticleType, true, true>::getNumFlopsPerKernelCall();
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
       << _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo) << endl;
  cout << "  Owned: " << _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned) << endl;
  cout << "  Halo : " << _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo) << endl;
  cout << "Standard Deviation of Homogeneity    : " << _homogeneity << endl;

  cout << fixed << setprecision(_floatStringPrecision);
  cout << "Measurements:" << endl;
  cout << timerToString("Time total      ", durationTotal, digitsTimeTotalNS, durationTotal);
  cout << timerToString("  Initialization", _timers.initialization.getTotalTime(), digitsTimeTotalNS, durationTotal);
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
  auto mfups = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned) * _iteration * 1e-6 /
               (_timers.forceUpdateTotal.getTotalTime() * 1e-9);  // 1e-9 for ns to s, 1e-6 for M in MFUP
  cout << "Tuning iterations: " << _numTuningIterations << " / " << _iteration << " = "
       << ((double)_numTuningIterations / _iteration * 100) << "%" << endl;
  cout << "MFUPs/sec    : " << mfups << endl;

  if (_configuration->dontMeasureFlops.value) {
    autopas::FlopCounterFunctor<ParticleType> flopCounterFunctor(_autoPasContainer->getCutoff());
    _autoPasContainer->iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * _iteration;
    // approximation for flops of verlet list generation
    if (_autoPasContainer->getContainerType() == autopas::ContainerOption::verletLists)
      flops += flopCounterFunctor.getDistanceCalculations() *
               decltype(flopCounterFunctor)::numFlopsPerDistanceCalculation *
               floor(_iteration / _configuration->verletRebuildFrequency.value);

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationSimulateSec << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }
}

