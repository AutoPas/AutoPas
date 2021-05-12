/**
 * @file MDFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexSingleNode.h"

#include "../BoundaryConditions.h"
#include "../Thermostat.h"
#include "../TimeDiscretization.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/utils/MemoryProfiler.h"

MDFlexSingleNode::MDFlexSingleNode(int dimensionCount, int argc, char **argv){
	MDFlexSimulation::initialize(dimensionCount, argc, argv);
}

void MDFlexSingleNode::run() {
  std::cout << std::endl << "Using " << autopas::autopas_get_max_threads() << " Threads" << std::endl;
  std::cout << "Starting simulation... " << std::endl;

  //this->_homogeneity = Simulation::calculateHomogeneity(autopas);
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
      calculatePositions();
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

    // invoke the force calculation with the functor specified in the configuration
    switch (_configuration->functorOption.value) {
      case MDFlexConfig::FunctorOption::lj12_6: {
        calculateForces<autopas::LJFunctor<ParticleType, _shifting, _mixing>>();
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_Globals: {
        calculateForces<autopas::LJFunctor<ParticleType, _shifting, _mixing, autopas::FunctorN3Modes::Both, true>>();
        break;
      }
      case MDFlexConfig::FunctorOption::lj12_6_AVX: {
        calculateForces<autopas::LJFunctorAVX<ParticleType, _shifting, _mixing>>();
        break;
      }
    }
    // only show memory usage in when the logger is set to debug
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << std::endl;
    }

    // only do time step related stuff when there actually is time-stepping
    if (_configuration->deltaT.value != 0) {
      _timers.velocityUpdate.start();
      TimeDiscretization::calculateVelocities(*_autoPasContainer, *_particlePropertiesLibrary, _configuration->deltaT.value);
      _timers.velocityUpdate.stop();

      // applying Velocity scaling with Thermostat:
      if (_configuration->useThermostat.value and (_iteration % _configuration->thermostatInterval.value) == 0) {
        _timers.thermostat.start();
        Thermostat::apply(*_autoPasContainer, *_particlePropertiesLibrary, _configuration->targetTemperature.value,
                          _configuration->deltaTemp.value);
        _timers.thermostat.stop();
      }
    }
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
    _configuration->initTemperature.value = Thermostat::calcTemperature(*_autoPasContainer, *_particlePropertiesLibrary);
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

void MDFlexSingleNode::initializeDomainDecomposition(int &dimensionCount){
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

template <class FunctorType>
void MDFlexSingleNode::calculateForces() {
  _timers.forceUpdateTotal.start();

  // pairwise forces
  _timers.forceUpdatePairwise.start();

  FunctorType functor{_autoPasContainer->getCutoff(), *_particlePropertiesLibrary};
  bool tuningIteration = _autoPasContainer->iteratePairwise(&functor);

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
  globalForces();
  _timers.forceUpdateGlobal.stop();
  _timers.forceUpdateTotal.stop();
}
