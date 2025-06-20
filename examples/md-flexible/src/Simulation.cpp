/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 01.03.2021
 */
#include "Simulation.h"

#include <algorithm>

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/SimilarityFunctions.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/utils/WrapOpenMP.h"

// Declare the main AutoPas class and the computeInteractions() methods with all used functors as extern template
// instantiation. They are instantiated in the respective cpp file inside the templateInstantiations folder.
//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeAutovec *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeAVX *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(LJFunctorTypeSVE *);
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC)
extern template bool autopas::AutoPas<ParticleType>::computeInteractions(ATFunctor *);
#endif
//! @endcond

#include <sys/ioctl.h>
#include <unistd.h>

#include <iostream>

#include "ParticleCommunicator.h"
#include "Thermostat.h"
#include "autopas/utils/MemoryProfiler.h"
#include "autopas/utils/WrapMPI.h"
#include "configuration/MDFlexConfig.h"

namespace {
/**
 * Tries to identify the width of the terminal where the simulation is running.
 * If no width can be identified, the function defaults to 80.
 * @return width of the terminal.
 */
size_t getTerminalWidth() {
  size_t terminalWidth = 0;
  // test all std pipes to get the current terminal width
  for (auto fd : {STDOUT_FILENO, STDIN_FILENO, STDERR_FILENO}) {
    if (isatty(fd)) {
      struct winsize w {};
      ioctl(fd, TIOCGWINSZ, &w);
      terminalWidth = w.ws_col;
      break;
    }
  }

  // if width is still zero try the environment variable COLUMNS
  if (terminalWidth == 0) {
    if (auto *terminalWidthCharArr = std::getenv("COLUMNS")) {
      // this pointer could be used to detect parsing errors via terminalWidthCharArr == end
      // but since we have a fallback further down we are ok if this fails silently.
      char *end{};
      terminalWidth = std::strtol(terminalWidthCharArr, &end, 10);
    }
  }

  // if all of the above fail fall back to a fixed width
  if (terminalWidth == 0) {
    // this seems to be the default width in most terminal windows
    terminalWidth = 80;
  }

  return terminalWidth;
}
}  // namespace

void Simulation::prepareAutopasContainer(std::shared_ptr<autopas::AutoPas<ParticleType>> &container, double cutoff,
                                         const std::string &outputSuffix) {
  container = std::make_shared<autopas::AutoPas<ParticleType>>(*_outputStream);
  container->setAllowedCellSizeFactors(*_configuration.cellSizeFactors.value);
  container->setAllowedContainers(_configuration.containerOptions.value);

  container->setAllowedInteractionTypeOptions(_configuration.getInteractionTypes());

  // Pairwise specific options
  container->setAllowedDataLayouts(_configuration.dataLayoutOptions.value, autopas::InteractionTypeOption::pairwise);
  container->setAllowedNewton3Options(_configuration.newton3Options.value, autopas::InteractionTypeOption::pairwise);
  container->setAllowedTraversals(_configuration.traversalOptions.value, autopas::InteractionTypeOption::pairwise);
  container->setAllowedLoadEstimators(_configuration.loadEstimatorOptions.value);
  // Triwise specific options
  container->setAllowedDataLayouts(_configuration.dataLayoutOptions3B.value, autopas::InteractionTypeOption::triwise);
  container->setAllowedNewton3Options(_configuration.newton3Options3B.value, autopas::InteractionTypeOption::triwise);
  container->setAllowedTraversals(_configuration.traversalOptions3B.value, autopas::InteractionTypeOption::triwise);
  // General options
  container->setBoxMin(_domainDecomposition->getLocalBoxMin());
  container->setBoxMax(_domainDecomposition->getLocalBoxMax());
  container->setCutoff(cutoff);
  container->setRelativeOptimumRange(_configuration.relativeOptimumRange.value);
  container->setMaxTuningPhasesWithoutTest(_configuration.maxTuningPhasesWithoutTest.value);
  container->setRelativeBlacklistRange(_configuration.relativeBlacklistRange.value);
  container->setEvidenceFirstPrediction(_configuration.evidenceFirstPrediction.value);
  container->setExtrapolationMethodOption(_configuration.extrapolationMethodOption.value);
  container->setNumSamples(_configuration.tuningSamples.value);
  container->setEarlyStoppingFactor(_configuration.earlyStoppingFactor.value);
  container->setMaxEvidence(_configuration.tuningMaxEvidence.value);
  container->setRuleFileName(_configuration.ruleFilename.value);
  container->setFuzzyRuleFileName(_configuration.fuzzyRuleFilename.value);
  container->setSelectorStrategy(_configuration.selectorStrategy.value);
  container->setTuningInterval(_configuration.tuningInterval.value);
  container->setTuningStrategyOption(_configuration.tuningStrategyOptions.value);
  container->setTuningMetricOption(_configuration.tuningMetricOption.value);
  container->setUseLOESSSmoothening(_configuration.useLOESSSmoothening.value);
  container->setEnergySensorOption(_configuration.energySensorOption.value);
  container->setMPITuningMaxDifferenceForBucket(_configuration.MPITuningMaxDifferenceForBucket.value);
  container->setMPITuningWeightForMaxDensity(_configuration.MPITuningWeightForMaxDensity.value);
  container->setVerletClusterSize(_configuration.verletClusterSize.value);
  container->setVerletRebuildFrequency(_configuration.verletRebuildFrequency.value);
  container->setVerletSkin(_configuration.verletSkinRadius.value);
  container->setAcquisitionFunction(_configuration.acquisitionFunctionOption.value);
  container->setUseTuningLogger(_configuration.useTuningLogger.value);
  container->setSortingThreshold(_configuration.sortingThreshold.value);
  container->setOutputSuffix(outputSuffix);

  container->init();
}

Simulation::Simulation(const MDFlexConfig &configuration,
                       std::shared_ptr<RegularGridDecomposition> &domainDecomposition)
    : _configuration(configuration),
      _domainDecomposition(domainDecomposition),
      _createVtkFiles(not configuration.vtkFileName.value.empty()),
      _vtkWriter(nullptr) {
  _timers.total.start();
  _timers.initialization.start();

  // only create the writer if necessary since this also creates the output dir
  if (_createVtkFiles) {
    _vtkWriter =
        std::make_shared<ParallelVtkWriter>(_configuration.vtkFileName.value, _configuration.vtkOutputFolder.value,
                                            std::to_string(_configuration.iterations.value).size());
  }

  const auto rank = _domainDecomposition->getDomainIndex();
  const auto *fillerBeforeSuffix =
      _configuration.outputSuffix.value.empty() or _configuration.outputSuffix.value.front() == '_' ? "" : "_";
  const auto *fillerAfterSuffix =
      _configuration.outputSuffix.value.empty() or _configuration.outputSuffix.value.back() == '_' ? "" : "_";
  const auto outputSuffix =
      "Rank" + std::to_string(rank) + fillerBeforeSuffix + _configuration.outputSuffix.value + fillerAfterSuffix;

  if (_configuration.logFileName.value.empty()) {
    _outputStream = &std::cout;
  } else {
    _logFile = std::make_shared<std::ofstream>();
    _logFile->open(_configuration.logFileName.value + "_" + outputSuffix);
    _outputStream = &(*_logFile);
  }

  if (_configuration.getInteractionTypes().empty()) {
    std::string functorName{};
    std::tie(_configuration.functorOption.value, functorName) =
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
        std::make_pair(MDFlexConfig::FunctorOption::lj12_6_AVX, "Lennard-Jones AVX Functor.");
#elif defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
        std::make_pair(MDFlexConfig::FunctorOption::lj12_6_SVE, "Lennard-Jones SVE Functor.");
#else
        std::make_pair(MDFlexConfig::FunctorOption::lj12_6, "Lennard-Jones AutoVec Functor.");
#endif
    _configuration.addInteractionType(autopas::InteractionTypeOption::pairwise);
    std::cout << "WARNING: No functor was specified. Defaulting to " << functorName << std::endl;
  }

  _useSecondAutopasInstanceForRespa = _configuration.useSecondAutpasInstance.value;

  const double cutoffToUse = _useSecondAutopasInstanceForRespa
                                 ? _configuration.cutoff.value
                                 : _configuration.cutoff.value * _configuration.cutoffFactorRespa.value;

  prepareAutopasContainer(_autoPasContainer, cutoffToUse, outputSuffix);

  if (_useSecondAutopasInstanceForRespa) {
    prepareAutopasContainer(_autoPasContainerRespa,
                            _configuration.cutoff.value * _configuration.cutoffFactorRespa.value, outputSuffix);
  } else {
    _domainDecomposition->setCutoff(cutoffToUse);
  }

  autopas::Logger::get()->set_level(_configuration.logLevel.value);

  // Throw an error if there is not more than one configuration to test in the search space but more than one tuning
  // phase is requested
  if (_autoPasContainer->searchSpaceIsTrivial() and _configuration.tuningPhases.value > 0) {
    throw std::runtime_error(
        "Search space must not be trivial if the simulation time is limited by the number tuning phases");
  }

  // Load particles from the config file (object generators, checkpoint)
  loadParticles();

  if (_configuration.useThermostat.value and _configuration.deltaT.value != 0) {
    if (_configuration.addBrownianMotion.value) {
      Thermostat::addBrownianMotion(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                    _configuration.initTemperature.value);
    }

    // Set the simulation directly to the desired initial temperature.
    Thermostat::apply(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                      _configuration.initTemperature.value, std::numeric_limits<double>::max());

    if (_configuration.respaStepSize.value > 1) {
      if (_configuration.thermostatInterval.value % _configuration.respaStepSize.value != 0) {
        throw std::runtime_error("The thermostat interval must be divisible by the respa-stepsize");
      }
    }
  }

  if (_configuration.respaStepSize.value > 1) {
    if (_configuration.iterations.value % _configuration.respaStepSize.value != 0) {
      throw std::runtime_error("The number of iterations must be divisible by the respa-stepsize");
    }
  }

  if (not _configuration.lutInputFile.value.empty() and _configuration.useApproxForceRespa.value) {
    throw std::runtime_error("Can not use LuT CG forces and approximate forces at the same time!");
  }

  if (_configuration.respaDistanceClassMode.value == autopas::DistanceClassOption::ibi and
      _configuration.lutInputFile.value.empty()) {
    throw std::runtime_error("Please specify a lut input file to be used with LuT functor in outer distance class!");
  }

  _timers.initialization.stop();

  // check if an ODF should be captured
  if (_configuration.odfNumBins.value > 0) {
    _odf = std::make_shared<ODF>(
        _autoPasContainer, 0, _configuration.odfRadius.value, _configuration.odfNumBins.value,
        _configuration.odfGuardArea.value,
        std::all_of(_configuration.boundaryOption.value.begin(), _configuration.boundaryOption.value.end(),
                    [](const auto &boundary) { return boundary == options::BoundaryTypeOption::periodic; })
            ? true
            : false);
  }

  // if an RDF should be created init it here
  const bool captureRDF = _configuration.rdfRadius.value > 0;
  if (captureRDF) {
    _rdfAA = std::make_shared<RDF>(
        _autoPasContainer, 0, _configuration.rdfRadius.value, _configuration.rdfNumBins.value,
        _configuration.rdfGuardArea.value,
        std::all_of(_configuration.boundaryOption.value.begin(), _configuration.boundaryOption.value.end(),
                    [](const auto &boundary) { return boundary == options::BoundaryTypeOption::periodic; })
            ? true
            : false);
    // if end iteration is not specified capture until end
    if (_configuration.rdfEndIteration.value == 0) {
      _configuration.rdfEndIteration.value = _configuration.iterations.value;
    }
  }

  // check if IBI should be used
  if (_configuration.ibiEquilibrateIterations.value > 0) {
    _rdfCG = std::make_shared<RDF>(
        _autoPasContainer, 0, _configuration.rdfRadius.value, _configuration.rdfNumBins.value,
        _configuration.rdfGuardArea.value,
        std::all_of(_configuration.boundaryOption.value.begin(), _configuration.boundaryOption.value.end(),
                    [](const auto &boundary) { return boundary == options::BoundaryTypeOption::periodic; })
            ? true
            : false);
  }

  if (not _configuration.lutInputFile.value.empty()) {
    std::cout << "Loading lookup table file: " << _configuration.lutInputFile.value << std::endl;
    _lut = std::make_shared<LookupTableType>();
    _lut->loadFromCSV(_configuration.lutInputFile.value);
    _lut->computeDerivatives();
    _cgSimulation = true;
  }

  if (_configuration.useApproxForceRespa.value) {
    _cgSimulation = true;
  }

  if (_configuration.respaDistanceClassMode.value != autopas::DistanceClassOption::disabled) {
    _distanceClassSimulation = true;
  }

  for (size_t i = 0; i < _configuration.getParticlePropertiesLibrary()->getNumberRegisteredSiteTypes(); ++i) {
    if (_configuration.getParticlePropertiesLibrary()->getCoulombEpsilon(i) > 0) {
      _applyCoulombFunctor = true;
    }
  }
}

void Simulation::finalize() {
  _timers.total.stop();
  autopas::AutoPas_MPI_Barrier(AUTOPAS_MPI_COMM_WORLD);

  logSimulationState();
  logMeasurements();
}

void Simulation::run() {
  _timers.simulate.start();

  bool equilibrate{false};
  size_t equilibrateIterations{0};
  size_t cgMeasureIterations{0};
  size_t ibiTrials{0};
  const size_t rdfNumBins = _configuration.rdfNumBins.value;
  bool ibiConvergenceReached{
      (_cgSimulation or _configuration.multiMultisiteModelsRespa.value or _distanceClassSimulation) ? true : false};
  bool firstRespaIterationSkipped{false};
  bool respaStarted{false};
  const bool ibiMeasureSimulation = _configuration.ibiEquilibrateIterations.value > 0;

  auto respaActive = _configuration.respaStepSize.value > 0;

  RotationalAnalysis rotationalAnalysis;
  if (_configuration.rotationalAnalysisLagSteps.value > 0) {
    rotationalAnalysis.setValues(_autoPasContainer, _configuration.rotationalAnalysisLagSteps.value,
                                 _configuration.rotationalAnalysisStepInterval.value);
  }

  while (needsMoreIterations()) {
    if (_createVtkFiles and _iteration % _configuration.vtkWriteFrequency.value == 0) {
      _timers.vtk.start();
      _vtkWriter->recordTimestep(_iteration, *_autoPasContainer, *_domainDecomposition);
      _timers.vtk.stop();
    }

    const bool isRespaIteration = respaActive ? _iteration % _configuration.respaStepSize.value == 0 : false;
    const bool nextIsRespaIteration = respaActive ? (_iteration + 1) % _configuration.respaStepSize.value == 0 : false;

    if ((_rdfCG and _iteration >= _configuration.rdfStartIteration.value and
         _iteration < _configuration.rdfEndIteration.value and
         _iteration % _configuration.rdfCaptureFreuency.value == 0) or
        (_rdfAA and not _rdfCG and _iteration >= _configuration.rdfStartIteration.value and
         _iteration < _configuration.rdfEndIteration.value and
         _iteration % _configuration.rdfCaptureFreuency.value == 0)) {
      _rdfAA->captureRDF();
    }

    // capture CG RDF
    if (_rdfCG and _iteration > _configuration.rdfEndIteration.value and not equilibrate and
        not ibiConvergenceReached) {
      _rdfCG->captureRDF();
      cgMeasureIterations++;
    }

    if (_odf and _iteration > _configuration.odfStartIteration.value and
        _iteration < _configuration.odfEndIteration.value and
        (_iteration % _configuration.odfCaptureFreuency.value == 0)) {
      _odf->captureODF();
    }

    _timers.computationalLoad.start();
    if (_configuration.deltaT.value != 0 and not _simulationIsPaused) {
      if (respaActive and isRespaIteration and ibiConvergenceReached) {
        // The first potential energy is not measured because it is not a full step.
        if (_distanceClassSimulation) {
          if (not respaStarted) {
            // do a first force evaluation
            if (_useSecondAutopasInstanceForRespa) {
              sendPositionsAndQuaternionsToRespaInstance();
            }
            if (_configuration.respaDistanceClassMode.value == autopas::DistanceClassOption::ibi) {
              updateInteractionForces(ForceType::IBIOuter);
            } else if (_configuration.respaDistanceClassMode.value == autopas::DistanceClassOption::fp) {
              updateInteractionForces(ForceType::FPOuter);
            } else if (_configuration.respaDistanceClassMode.value == autopas::DistanceClassOption::cgmol) {
              updateInteractionForces(ForceType::CGMolOuter);
            }
            if (_useSecondAutopasInstanceForRespa) {
              sendBackForcesAndTorquesFromRespaInstance();
            }
          }
          updateVelocities(/*resetForce*/ true, RespaIterationType::OuterStep);
#if MD_FLEXIBLE_MODE == MULTISITE
          if (_configuration.multiMultisiteModelsRespa.value) {
            updateAngularVelocities(/*resetTorques*/ true, RespaIterationType::OuterStep,
                                    _configuration.respaMoleculeTypes.value[0],
                                    _configuration.respaMoleculeTypes.value[1]);
          } else {
            updateAngularVelocities(/*resetTorques*/ true, RespaIterationType::OuterStep, 0, 0);
          }
#endif
          respaStarted = true;
        } else {
          if (firstRespaIterationSkipped) {
            if (not respaStarted) {
              // do a first force evaluation
              updateInteractionForces(ForceType::FullParticle);
              if (not _configuration.multiMultisiteModelsRespa.value) {
                updateInteractionForces(ForceType::CoarseGrain, true);
              }
            }
            updateVelocities(/*resetForce*/ true, RespaIterationType::OuterStep);
#if MD_FLEXIBLE_MODE == MULTISITE
            if (_configuration.multiMultisiteModelsRespa.value) {
              updateAngularVelocities(/*resetTorques*/ true, RespaIterationType::OuterStep,
                                      _configuration.respaMoleculeTypes.value[0],
                                      _configuration.respaMoleculeTypes.value[1]);

            } else {
              updateAngularVelocities(/*resetTorques*/ false, RespaIterationType::OuterStep, 0, 0);
            }
#endif

            respaStarted = true;
          } else {
            firstRespaIterationSkipped = true;
          }
        }
      }

      updatePositionsAndResetForces();
#if MD_FLEXIBLE_MODE == MULTISITE
      // updateQuaternions(/*resetTorques*/ not(respaStarted and not nextIsRespaIteration));
      updateQuaternions();
#endif

      _timers.updateContainer.start();
      auto emigrants = _autoPasContainer->updateContainer();
      _timers.updateContainer.stop();

      const auto computationalLoad = static_cast<double>(_timers.computationalLoad.stop());

      // periodically resize box for MPI load balancing
      if (_iteration % _configuration.loadBalancingInterval.value == 0) {
        _timers.loadBalancing.start();
        _domainDecomposition->update(computationalLoad);
        auto additionalEmigrants = _autoPasContainer->resizeBox(_domainDecomposition->getLocalBoxMin(),
                                                                _domainDecomposition->getLocalBoxMax());
        // If the boundaries shifted, particles that were thrown out by updateContainer() previously might now be in the
        // container again.
        // Reinsert emigrants if they are now inside the domain and mark local copies as dummy,
        // so that remove_if can erase them after.
        const auto &boxMin = _autoPasContainer->getBoxMin();
        const auto &boxMax = _autoPasContainer->getBoxMax();
        _autoPasContainer->addParticlesIf(emigrants, [&](auto &p) {
          if (autopas::utils::inBox(p.getR(), boxMin, boxMax)) {
            // This only changes the ownership state in the emigrants vector, not in AutoPas
            p.setOwnershipState(autopas::OwnershipState::dummy);
            return true;
          }
          return false;
        });

        emigrants.erase(std::remove_if(emigrants.begin(), emigrants.end(), [&](const auto &p) { return p.isDummy(); }),
                        emigrants.end());

        emigrants.insert(emigrants.end(), additionalEmigrants.begin(), additionalEmigrants.end());
        _timers.loadBalancing.stop();
      }

      _timers.migratingParticleExchange.start();
      _domainDecomposition->exchangeMigratingParticles(*_autoPasContainer, emigrants);
      _timers.migratingParticleExchange.stop();

      _timers.reflectParticlesAtBoundaries.start();
      _domainDecomposition->reflectParticlesAtBoundaries(*_autoPasContainer,
                                                         *_configuration.getParticlePropertiesLibrary());
      _timers.reflectParticlesAtBoundaries.stop();

      _timers.haloParticleExchange.start();
      _domainDecomposition->exchangeHaloParticles(
          *_autoPasContainer, (_useSecondAutopasInstanceForRespa and _distanceClassSimulation)
                                  ? _configuration.cutoff.value
                                  : _configuration.cutoff.value * _configuration.cutoffFactorRespa.value);
      _timers.haloParticleExchange.stop();

      _timers.computationalLoad.start();
    }

    if (_distanceClassSimulation and respaActive) {
      const auto potentialEnergyAndVirial = updateInteractionForces(ForceType::FPInner);
      if (nextIsRespaIteration and mdFlexibleTypeDefs::calcGlobals) {
        _potentialEnergyInnerLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).first);
        _virialInnerLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).second);
      }
    } else {
      if ((not respaActive and not _cgSimulation and not ibiMeasureSimulation) or
          (_iteration <= _configuration.rdfEndIteration.value and not _cgSimulation and not respaActive) or
          (respaActive and _configuration.multiMultisiteModelsRespa.value)) {
        const auto potentialEnergyAndVirial = updateInteractionForces(ForceType::FullParticle);
        if (((respaActive and nextIsRespaIteration) or (not respaActive)) and mdFlexibleTypeDefs::calcGlobals) {
          _potentialEnergyInnerLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).first);
          _virialInnerLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).second);
        }
      } else {
        const auto potentialEnergyAndVirial = updateInteractionForces(ForceType::CoarseGrain);
        if (((respaActive and nextIsRespaIteration and respaStarted) or
             (_cgSimulation and not respaActive) and mdFlexibleTypeDefs::calcGlobals)) {
          _potentialEnergyInnerLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).first);
          _virialInnerLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).second);
        }
      }
    }

    if (_configuration.pauseSimulationDuringTuning.value) {
      // If PauseSimulationDuringTuning is enabled we need to update the _simulationIsPaused flag
      updateSimulationPauseState();
    }

    if (_configuration.deltaT.value != 0 and not _simulationIsPaused) {
      updateVelocities(/*resetForces*/ respaActive and nextIsRespaIteration and ibiConvergenceReached);
#if MD_FLEXIBLE_MODE == MULTISITE
      if (_configuration.multiMultisiteModelsRespa.value) {
        updateAngularVelocities(/*resetTorques*/ respaActive and nextIsRespaIteration and ibiConvergenceReached,
                                RespaIterationType::InnerStep, _configuration.respaMoleculeTypes.value[0],
                                _configuration.respaMoleculeTypes.value[1]);
      } else if (_distanceClassSimulation and not _configuration.multiMultisiteModelsRespa.value) {
        updateAngularVelocities(/*resetTorques*/ respaActive and nextIsRespaIteration and ibiConvergenceReached,
                                RespaIterationType::InnerStep, 0, 0);
      } else {
        updateAngularVelocities(false, RespaIterationType::InnerStep, 0, 0);
      }

#endif
      if (_configuration.useThermostat.value) {
        if ((not respaActive or not respaStarted) and
            ((_iteration + 1) % _configuration.thermostatInterval.value == 0)) {
          updateThermostat(true);
        }
      }

      if (not respaActive) {
        if (_potentialEnergyInnerLoop.size() > 0) {
          _kineticEnergy.push_back(calculateKineticEnergy());
          _totalEnergy.push_back(_potentialEnergyInnerLoop.back() + _kineticEnergy.back());
        }
      }

      if (respaActive and nextIsRespaIteration and ibiConvergenceReached and respaStarted) {
        std::optional<std::pair<double, double>> potentialEnergyAndVirial;
        if (_distanceClassSimulation) {
          if (_useSecondAutopasInstanceForRespa) {
            sendPositionsAndQuaternionsToRespaInstance();
          }
          if (_configuration.respaDistanceClassMode.value == autopas::DistanceClassOption::ibi) {
            potentialEnergyAndVirial = updateInteractionForces(ForceType::IBIOuter);
          } else if (_configuration.respaDistanceClassMode.value == autopas::DistanceClassOption::fp) {
            potentialEnergyAndVirial = updateInteractionForces(ForceType::FPOuter);
          } else if (_configuration.respaDistanceClassMode.value == autopas::DistanceClassOption::cgmol) {
            potentialEnergyAndVirial = updateInteractionForces(ForceType::CGMolOuter);
          }
          if (_useSecondAutopasInstanceForRespa) {
            sendBackForcesAndTorquesFromRespaInstance();
          }
          updateVelocities(false, RespaIterationType::OuterStep);
#if MD_FLEXIBLE_MODE == MULTISITE
          if (_configuration.multiMultisiteModelsRespa.value) {
            updateAngularVelocities(false, RespaIterationType::OuterStep, _configuration.respaMoleculeTypes.value[0],
                                    _configuration.respaMoleculeTypes.value[1]);
          } else {
            updateAngularVelocities(false, RespaIterationType::OuterStep, 0, 0);
          }
#endif
        } else {
          potentialEnergyAndVirial = updateInteractionForces(ForceType::FullParticle);
          if (not _configuration.multiMultisiteModelsRespa.value) {
            potentialEnergyAndVirial = updateInteractionForces(ForceType::CoarseGrain, true);
          }
          updateVelocities(false, RespaIterationType::OuterStep);
#if MD_FLEXIBLE_MODE == MULTISITE
          if (_configuration.multiMultisiteModelsRespa.value) {
            updateAngularVelocities(false, RespaIterationType::OuterStep, _configuration.respaMoleculeTypes.value[0],
                                    _configuration.respaMoleculeTypes.value[1]);
          } else {
            updateAngularVelocities(false, RespaIterationType::OuterStep, 0, 0);
          }

#endif
        }

        if (mdFlexibleTypeDefs::calcGlobals) {
          _potentialEnergyOuterLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).first);
          _virialOuterLoop.push_back(potentialEnergyAndVirial.value_or(std::pair(0, 0)).second);
        }

        if (_configuration.useThermostat.value) {
          if ((_iteration + 1) % _configuration.thermostatInterval.value == 0) {
            updateThermostat(/*skipIterationCheck*/ true);
          }
        }

        if (_potentialEnergyInnerLoop.size() > 0 and _potentialEnergyOuterLoop.size() > 0) {
          _kineticEnergy.push_back(calculateKineticEnergy());
          _totalEnergy.push_back(_potentialEnergyInnerLoop.back() + _potentialEnergyOuterLoop.back() +
                                 _kineticEnergy.back());
        }
      }
    }
    _timers.computationalLoad.stop();

    if (_configuration.rotationalAnalysisLagSteps.value > 0 and
        _iteration >= _configuration.rotationalAnalysisStartIteration.value and
        _iteration <= _configuration.rotationalAnalysisEndIteration.value) {
      rotationalAnalysis.recordOrientations(_iteration);
    }

    if (not _simulationIsPaused) {
      ++_iteration;
    }

    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::cout << "Current Memory usage on rank " << _domainDecomposition->getDomainIndex() << ": "
                << autopas::memoryProfiler::currentMemoryUsage() << " kB\n";
    }

    if (_domainDecomposition->getDomainIndex() == 0) {
      auto [maxIterationsEstimate, maxIterationsIsPrecise] = estimateNumberOfIterations();
      if (not _configuration.dontShowProgressBar.value) {
        printProgress(_iteration, maxIterationsEstimate, maxIterationsIsPrecise);
      }
    }

    // capture the last RDF from original simulation
    if (_rdfCG and _iteration == _configuration.rdfEndIteration.value) {
      _rdfAA->captureRDF();
      _rdfAA->computeFinalRDF();

      _rdfAA->writeToCSV(_configuration.rdfOutputFolder.value,
                         _configuration.rdfFileName.value + "_" + std::to_string(_iteration));

      // generate initial cg_potential
      _lut = std::make_shared<LookupTableType>();

      // this assumes normalized values for the temperature. Therefore the value equals temperature * boltzmannConstant
      const double kT = _configuration.targetTemperature.value;

      std::vector<std::pair<double, double>> initialPotential(rdfNumBins, std::pair<double, double>(0, 0));
      int minIdx = -1;
      for (int i = 0; i < rdfNumBins; i++) {
        initialPotential[i].first = std::pow(
            _configuration.rdfRadius.value / static_cast<double>(rdfNumBins) * static_cast<double>(i + 1), 2.0);
        if (_rdfAA->getFinalRDF()[i].second > 0.0) {
          initialPotential[i].second = -kT * std::log(_rdfAA->getFinalRDF()[i].second);
        } else {
          minIdx = i;
        }
      }
      if (minIdx != -1) {
        const double dy = initialPotential[minIdx + 2].second - initialPotential[minIdx + 1].second;
        double yVal = initialPotential[minIdx + 1].second;
        for (int idx = minIdx + 1; idx >= 0; idx--) {
          initialPotential[idx].second = yVal;
          yVal -= dy;
        }
      }

      _lut->updateTable(initialPotential);
      _lut->computeDerivatives();
      _lut->writeToCSV(_configuration.lutOutputFolder.value, _configuration.lutFileName.value + "_0");

      equilibrate = true;
    }

    if (_odf and _iteration == _configuration.odfEndIteration.value) {
      _odf->captureODF();
      _odf->computeFinalODF();
      _odf->writeToCSV(_configuration.odfOutputFolder.value,
                       _configuration.odfFileName.value + "_" + std::to_string(_iteration));
    }

    if (_configuration.rotationalAnalysisLagSteps.value > 0 and
        _iteration == _configuration.rotationalAnalysisEndIteration.value) {
      rotationalAnalysis.writeResults(_configuration.rotationalAnalysisOutputFolder.value,
                                      _configuration.rotationalAnalysisFilename.value);
    }

    if (_rdfAA and not _rdfCG and _iteration == _configuration.rdfEndIteration.value) {
      _rdfAA->captureRDF();
      _rdfAA->computeFinalRDF();

      _rdfAA->writeToCSV(_configuration.rdfOutputFolder.value,
                         _configuration.rdfFileName.value + "_" + std::to_string(_iteration));
    }

    if (_rdfCG and _iteration > _configuration.rdfEndIteration.value and equilibrate) {
      equilibrateIterations++;
      if (equilibrateIterations == _configuration.ibiEquilibrateIterations.value) {
        equilibrate = false;
        equilibrateIterations = 0;
      }
    }

    // get the CG RDF and compare it to the AA RDF
    if (_rdfCG and
        cgMeasureIterations == (_configuration.rdfEndIteration.value - _configuration.rdfStartIteration.value) and
        not ibiConvergenceReached) {
      _rdfCG->captureRDF();
      _rdfCG->computeFinalRDF();

      _rdfCG->writeToCSV(_configuration.rdfOutputFolder.value,
                         _configuration.rdfFileName.value + "_" + std::to_string(_iteration));

      const auto rdfCGValues = _rdfCG->getFinalRDF();
      const auto rdfRefValues = _rdfAA->getFinalRDF();

      // this assumes normalized values for the temperature. Therefore the value equals temperature * boltzmannConstant
      const auto kT = _configuration.targetTemperature.value;

      const auto tolerance = _configuration.ibiConvergenceThreshold.value;

      const auto updateAlpha = _configuration.ibiUpdateAlpha.value;

      // Update the CG potential
      int minIdx = -1;
      for (int i = 0; i < rdfNumBins; i++) {
        if (rdfRefValues[i].second > 0.0 and rdfCGValues[i].second > 0.0) {
          double update = updateAlpha * kT * std::log(rdfCGValues[i].second / rdfRefValues[i].second);
          (*_lut)[i] += update;
        } else {
          minIdx = i;
        }
      }
      // add correction for values where RDF is 0
      if (minIdx != -1) {
        const double dy = (*_lut)[minIdx + 2] - (*_lut)[minIdx + 1];
        double yVal = (*_lut)[minIdx + 1];
        for (int idx = minIdx + 1; idx >= 0; idx--) {
          (*_lut)[idx] = yVal;
          yVal -= dy;
        }
      }

      // convergence check
      double sum0{0};
      double sum1{0};
      for (size_t i = 0; i < rdfNumBins; ++i) {
        sum0 += std::abs(rdfCGValues[i].second - rdfRefValues[i].second);
        sum1 += (std::abs(rdfCGValues[i].second) + std::abs(rdfRefValues[i].second));
      }
      double convergenceResult = 1.0 - sum0 / sum1;

      ibiTrials++;

      _lut->writeToCSV(_configuration.lutOutputFolder.value,
                       _configuration.lutFileName.value + "_" + std::to_string(ibiTrials));

      _lut->computeDerivatives();

      // continue with CG fitting or stop if convergence reached
      if (convergenceResult < tolerance) {
        _rdfCG->reset();
        cgMeasureIterations = 0;
        equilibrate = true;
        std::cout << "IBI: still fitting! convergenceResult: " << convergenceResult << std::endl;
      } else {
        ibiConvergenceReached = true;
        std::cout << "IBI: convergence reached!" << std::endl;
        exit(0);
      }
    }
  }
  _timers.simulate.stop();

  // Record last state of simulation.
  if (_createVtkFiles) {
    _vtkWriter->recordTimestep(_iteration, *_autoPasContainer, *_domainDecomposition);
  }

  if (mdFlexibleTypeDefs::calcGlobals) {
    bool allVectorSizesEqual = false;
    std::array<size_t, 6> sizesOfVectors = {
        _potentialEnergyOuterLoop.size(), _potentialEnergyInnerLoop.size(), _kineticEnergy.size(), _totalEnergy.size(),
        _virialInnerLoop.size(),          _virialOuterLoop.size()};

    if (not respaActive) {
      sizesOfVectors = {_potentialEnergyInnerLoop.size(),
                        _potentialEnergyInnerLoop.size(),
                        _kineticEnergy.size(),
                        _totalEnergy.size(),
                        _virialInnerLoop.size(),
                        _virialInnerLoop.size()};
    }

    allVectorSizesEqual = std::all_of(sizesOfVectors.begin(), sizesOfVectors.end(),
                                      [&sizesOfVectors](const auto value) { return value == sizesOfVectors[0]; });
    if (!allVectorSizesEqual) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "_potentialEnergyInnerLoop.size() = " + std::to_string(_potentialEnergyInnerLoop.size()) +
          " _kineticEnergy.size() = " + std::to_string(_kineticEnergy.size()) + " _totalEnergy.size() = " +
          std::to_string(_totalEnergy.size()) + " : somethong is wrong with the globals calculation");
    }

    // prepare JSON-output
    std::ostringstream json;
    json << std::setprecision(15);
    json << "{\n";

    auto writeVector = [&json](const std::string &name, const std::vector<double> &vec) {
      json << "  \"" << name << "\": [";
      for (size_t i = 0; i < vec.size(); ++i) {
        json << vec[i];
        if (i != vec.size() - 1) {
          json << ", ";
        }
      }
      json << "],\n";
    };

    // print to stdout
    std::cout << "potential energy inner loop: [";
    for (size_t i = 0; i < _potentialEnergyInnerLoop.size(); ++i) {
      std::cout << std::setprecision(15) << _potentialEnergyInnerLoop[i];
      if (i != _potentialEnergyInnerLoop.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;

    writeVector("potentialEnergyInnerLoop", _potentialEnergyInnerLoop);

    std::cout << "virial inner loop: [";
    for (size_t i = 0; i < _virialInnerLoop.size(); ++i) {
      std::cout << std::setprecision(15) << _virialInnerLoop[i];
      if (i != _virialInnerLoop.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;

    writeVector("virialInnerLoop", _virialInnerLoop);

    if (_potentialEnergyOuterLoop.size() > 0) {
      std::cout << "potential energy outer loop: [";
      for (size_t i = 0; i < _potentialEnergyOuterLoop.size(); ++i) {
        std::cout << std::setprecision(15) << _potentialEnergyOuterLoop[i];
        if (i != _potentialEnergyOuterLoop.size() - 1) {
          std::cout << ", ";
        }
      }
      std::cout << "]" << std::endl;

      writeVector("potentialEnergyOuterLoop", _potentialEnergyOuterLoop);
    }

    if (_virialOuterLoop.size() > 0) {
      std::cout << "virial outer loop: [";
      for (size_t i = 0; i < _virialOuterLoop.size(); ++i) {
        std::cout << std::setprecision(15) << _virialOuterLoop[i];
        if (i != _virialOuterLoop.size() - 1) {
          std::cout << ", ";
        }
      }
      std::cout << "]" << std::endl;

      writeVector("virialOuterLoop", _virialOuterLoop);
    }

    std::vector<double> combinedPotentialEnergy;
    if (_potentialEnergyOuterLoop.size() == _potentialEnergyInnerLoop.size()) {
      std::cout << "potential energy combined: [";
      for (size_t i = 0; i < _potentialEnergyInnerLoop.size(); ++i) {
        const auto sum = _potentialEnergyInnerLoop[i] + _potentialEnergyOuterLoop[i];
        combinedPotentialEnergy.push_back(sum);
        std::cout << std::setprecision(15) << sum;
        if (i != _potentialEnergyInnerLoop.size() - 1) {
          std::cout << ", ";
        }
      }
      std::cout << "]" << std::endl;

      writeVector("potentialEnergyCombined", combinedPotentialEnergy);
    }

    std::vector<double> combinedVirial;
    if (_virialOuterLoop.size() == _virialInnerLoop.size()) {
      std::cout << "virial combined: [";
      for (size_t i = 0; i < _virialInnerLoop.size(); ++i) {
        const auto sum = _virialInnerLoop[i] + _virialOuterLoop[i];
        combinedVirial.push_back(sum);
        std::cout << std::setprecision(15) << sum;
        if (i != _virialInnerLoop.size() - 1) {
          std::cout << ", ";
        }
      }
      std::cout << "]" << std::endl;

      writeVector("combinedVirial", combinedVirial);
    }

    std::cout << "kinetic energy: [";
    for (size_t i = 0; i < _kineticEnergy.size(); ++i) {
      std::cout << std::setprecision(15) << _kineticEnergy[i];
      if (i != _kineticEnergy.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;

    std::cout << "total energy: [";
    for (size_t i = 0; i < _totalEnergy.size(); ++i) {
      std::cout << std::setprecision(15) << _totalEnergy[i];
      if (i != _totalEnergy.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;

    writeVector("kineticEnergy", _kineticEnergy);
    writeVector("totalEnergy", _totalEnergy);

    const auto numSteps = _totalEnergy.size();
    const auto avgKineticEnergy =
        std::reduce(_kineticEnergy.begin(), _kineticEnergy.end()) / static_cast<double>(numSteps);
    const auto avgTotalEnergy = std::reduce(_totalEnergy.begin(), _totalEnergy.end()) / static_cast<double>(numSteps);

    double sum = 0.0;
    for (const auto eT : _totalEnergy) {
      sum += std::abs(eT - avgTotalEnergy);
    }

    double rvite = sum / (avgKineticEnergy * static_cast<double>(_configuration.iterations.value));

    std::cout << "RelativeVariationInTrueEnergy: " << std::setprecision(15) << rvite << std::endl;
    json << "  \"RelativeVariationInTrueEnergy\": " << rvite << "\n";
    json << "}\n";

    // Write JSON
    std::filesystem::create_directories(_configuration.statisticsOutputFolder.value);
    std::ofstream outFile(_configuration.statisticsOutputFolder.value + "/" +
                          _configuration.statisticsOutputFilename.value + ".json");
    if (outFile.is_open()) {
      outFile << json.str();
      outFile.close();
    } else {
      std::cerr << "Failed to open statistics.json for writing!" << std::endl;
    }
  }
}

void Simulation::sendPositionsAndQuaternionsToRespaInstance() {
  _timers.secondAutopasInstanceSync.start();
  if (_autoPasContainerRespa->getNumberOfParticles() != _autoPasContainer->getNumberOfParticles()) {
    // do a first sync at the beginning of the simulation (copy particles to second autopas instance)
    _autoPasContainerRespa->deleteAllParticles();
    for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      _autoPasContainerRespa->addParticle(*particle);
    }
  } else {
    // 1. update positions and quaternions and reset forces (aka sync positions, quaternions from
    // inner-respa-loop-instance and set the forces to 0 or the global force)
    std::unordered_map<size_t, std::pair<std::array<double, 3UL>, std::array<double, 4UL>>>
        positionsAndQuaternionsOwned;
    positionsAndQuaternionsOwned.reserve(_autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned));

    for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      positionsAndQuaternionsOwned.emplace(particle->getID(), std::pair(particle->getR(), particle->getQuaternion()));
    }

    for (auto particle = _autoPasContainerRespa->begin(autopas::IteratorBehavior::owned); particle.isValid();
         ++particle) {
      const auto pR = positionsAndQuaternionsOwned[particle->getID()];
      particle->setR(pR.first);
      particle->setQuaternion(pR.second);
      particle->setF({0, 0, 0});
      particle->setTorque({0, 0, 0});
    }
  }

  // update the container (removes halos)
  const auto leavingParticles = _autoPasContainerRespa->updateContainer();

  // adjust box size to be the same as the other autopas instance
  auto additionalEmigrants =
      _autoPasContainerRespa->resizeBox(_domainDecomposition->getLocalBoxMin(), _domainDecomposition->getLocalBoxMax());

  // sanity check
  if (leavingParticles.size() != 0 or additionalEmigrants.size() != 0) {
    throw autopas::utils::ExceptionHandler::AutoPasException(
        "leavingParticles or additionalEmigrants in _autoPasContainerRespa is not 0 after particle sync and "
        "updateContainer");
  }

  // exchange and reflect happend already in the primary autopas instance

  // recreate halos with larger cutoff
  _domainDecomposition->exchangeHaloParticles(*_autoPasContainerRespa,
                                              _configuration.cutoff.value * _configuration.cutoffFactorRespa.value);
  _timers.secondAutopasInstanceSync.stop();
}

void Simulation::sendBackForcesAndTorquesFromRespaInstance() {
  _timers.secondAutopasInstanceSync.start();
  std::unordered_map<size_t, std::pair<std::array<double, 3UL>, std::array<double, 3UL>>> forcesAndTorquesOwned;
  forcesAndTorquesOwned.reserve(_autoPasContainerRespa->getNumberOfParticles(autopas::IteratorBehavior::owned));

  for (auto particle = _autoPasContainerRespa->begin(autopas::IteratorBehavior::owned); particle.isValid();
       ++particle) {
    forcesAndTorquesOwned.emplace(particle->getID(), std::pair(particle->getF(), particle->getTorque()));
  }

  for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto pR = forcesAndTorquesOwned[particle->getID()];
    particle->setF(pR.first);
    particle->setTorque(pR.second);
  }
  _timers.secondAutopasInstanceSync.stop();
}

std::tuple<size_t, bool> Simulation::estimateNumberOfIterations() const {
  if (_configuration.tuningPhases.value > 0) {
    const size_t configsTestedPerTuningPhase = [&]() {
      if (std::any_of(_configuration.tuningStrategyOptions.value.begin(),
                      _configuration.tuningStrategyOptions.value.end(), [](const auto &stratOpt) {
                        return stratOpt == autopas::TuningStrategyOption::bayesianSearch or
                               stratOpt == autopas::TuningStrategyOption::bayesianClusterSearch or
                               stratOpt == autopas::TuningStrategyOption::randomSearch;
                      })) {
        return static_cast<size_t>(_configuration.tuningMaxEvidence.value);
      } else {
        // @TODO: this can be improved by considering the tuning strategy
        //      or averaging number of iterations per tuning phase and dynamically adapt prediction

        // This estimate is only valid for full search and no restrictions on the cartesian product.
        // add static to only evaluate this once
        const size_t searchSpaceSizePairwise =
            _configuration.getInteractionTypes().count(autopas::InteractionTypeOption::pairwise) == 0
                ? 0
                : autopas::SearchSpaceGenerators::cartesianProduct(
                      _configuration.containerOptions.value, _configuration.traversalOptions.value,
                      _configuration.loadEstimatorOptions.value, _configuration.dataLayoutOptions.value,
                      _configuration.newton3Options.value, _configuration.cellSizeFactors.value.get(),
                      autopas::InteractionTypeOption::pairwise)
                      .size();

        const size_t searchSpaceSizeTriwise =
            _configuration.getInteractionTypes().count(autopas::InteractionTypeOption::triwise) == 0
                ? 0
                : autopas::SearchSpaceGenerators::cartesianProduct(
                      _configuration.containerOptions.value, _configuration.traversalOptions3B.value,
                      _configuration.loadEstimatorOptions.value, _configuration.dataLayoutOptions3B.value,
                      _configuration.newton3Options3B.value, _configuration.cellSizeFactors.value.get(),
                      autopas::InteractionTypeOption::triwise)
                      .size();

        return std::max(searchSpaceSizePairwise, searchSpaceSizeTriwise);
      }
    }();
    // non-tuning iterations + tuning iterations + one iteration after last phase
    auto estimate =
        (_configuration.tuningPhases.value - 1) * _configuration.tuningInterval.value +
        (_configuration.tuningPhases.value * _configuration.tuningSamples.value * configsTestedPerTuningPhase) + 1;
    return {estimate, false};
  } else {
    return {_configuration.iterations.value, true};
  }
}

void Simulation::printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise) {
  // percentage of iterations complete
  const double fractionDone = static_cast<double>(iterationProgress) / static_cast<double>(maxIterations);

  // length of the number of maxIterations
  const auto numCharsOfMaxIterations = static_cast<int>(std::to_string(maxIterations).size());

  // trailing information string
  std::stringstream info;
  info << std::setw(3) << std::round(fractionDone * 100) << "% " << std::setw(numCharsOfMaxIterations)
       << iterationProgress << "/";
  if (not maxIsPrecise) {
    info << "~";
  }
  info << maxIterations;

  // actual progress bar
  std::stringstream progressbar;
  progressbar << "[";
  // get current terminal width
  const auto terminalWidth = getTerminalWidth();

  // the bar should fill the terminal window so subtract everything else (-2 for "] ")
  const int maxBarWidth = static_cast<int>(terminalWidth - info.str().size() - progressbar.str().size() - 2ul);
  // sanity check for underflow
  if (maxBarWidth > terminalWidth) {
    std::cerr << "Warning! Terminal width appears to be too small or could not be read. Disabling progress bar."
              << std::endl;
    _configuration.dontShowProgressBar.value = true;
    return;
  }
  const auto barWidth =
      std::max(std::min(static_cast<decltype(maxBarWidth)>(maxBarWidth * (fractionDone)), maxBarWidth), 1);
  // don't print arrow tip if >= 100%
  if (iterationProgress >= maxIterations) {
    progressbar << std::string(barWidth, '=');
  } else {
    progressbar << std::string(barWidth - 1, '=') << '>' << std::string(maxBarWidth - barWidth, ' ');
  }
  progressbar << "] ";
  // clear current line (=delete previous progress bar)
  std::cout << std::string(terminalWidth, '\r');
  // print everything
  std::cout << progressbar.str() << info.str() << std::flush;
}

std::string Simulation::timerToString(const std::string &name, long timeNS, int numberWidth, long maxTime) {
  // only print timers that were actually used
  if (timeNS == 0) {
    return "";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(_floatStringPrecision) << name << " : " << std::setw(numberWidth) << std::right
     << timeNS
     << " ns ("
     // min width of the representation of seconds is numberWidth - 9 (from conversion) + 4 (for dot and digits after)
     << std::setw(numberWidth - 5) << ((double)timeNS * 1e-9) << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeNS / (double)maxTime * 100) << "%";
  }
  ss << "\n";
  return ss.str();
}

void Simulation::updatePositionsAndResetForces() {
  _timers.positionUpdate.start();
  TimeDiscretization::calculatePositionsAndResetForces(
      *_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()), _configuration.deltaT.value,
      _configuration.globalForce.value, _configuration.fastParticlesThrow.value);
  _timers.positionUpdate.stop();
}

void Simulation::updateQuaternions(bool resetTorques) {
  _timers.quaternionUpdate.start();
  TimeDiscretization::calculateQuaternionsAndResetTorques(
      *_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()), _configuration.deltaT.value,
      _configuration.globalForce.value, resetTorques);
  _timers.quaternionUpdate.stop();
}

std::optional<std::pair<double, double>> Simulation::updateInteractionForces(ForceType forceTypeToCalculate,
                                                                             bool subtractForces) {
  _timers.forceUpdateTotal.start();

  _previousIterationWasTuningIteration = _currentIterationIsTuningIteration;
  _currentIterationIsTuningIteration = false;
  long timeIteration = 0;
  double potentialEnergy = 0;
  double virial = 0;

  // Calculate pairwise forces
  if (_configuration.getInteractionTypes().count(autopas::InteractionTypeOption::pairwise)) {
    _timers.forceUpdatePairwise.start();
    std::tuple<bool, std::optional<double>, std::optional<double>> wasTuningIterationPotentialEnergyVirial =
        calculatePairwiseForces(forceTypeToCalculate, subtractForces);
    _currentIterationIsTuningIteration =
        (_currentIterationIsTuningIteration | std::get<0>(wasTuningIterationPotentialEnergyVirial));
    if (std::get<1>(wasTuningIterationPotentialEnergyVirial).has_value()) {
      potentialEnergy += std::get<1>(wasTuningIterationPotentialEnergyVirial).value();
    }
    if (std::get<2>(wasTuningIterationPotentialEnergyVirial).has_value()) {
      virial += std::get<2>(wasTuningIterationPotentialEnergyVirial).value();
    }
    timeIteration += _timers.forceUpdatePairwise.stop();
  }
  // Calculate triwise forces
  if (_configuration.getInteractionTypes().count(autopas::InteractionTypeOption::triwise)) {
    _timers.forceUpdateTriwise.start();
    std::tuple<bool, std::optional<double>, std::optional<double>> wasTuningIterationPotentialEnergyVirial =
        calculateTriwiseForces(forceTypeToCalculate);
    _currentIterationIsTuningIteration =
        (_currentIterationIsTuningIteration | std::get<0>(wasTuningIterationPotentialEnergyVirial));
    if (std::get<1>(wasTuningIterationPotentialEnergyVirial).has_value()) {
      potentialEnergy += std::get<1>(wasTuningIterationPotentialEnergyVirial).value();
    }
    if (std::get<2>(wasTuningIterationPotentialEnergyVirial).has_value()) {
      virial += std::get<2>(wasTuningIterationPotentialEnergyVirial).value();
    }
    timeIteration += _timers.forceUpdateTriwise.stop();
  }

  // count time spent for tuning
  if (_currentIterationIsTuningIteration) {
    _timers.forceUpdateTuning.addTime(timeIteration);
    ++_numTuningIterations;
  } else {
    _timers.forceUpdateNonTuning.addTime(timeIteration);
    // if the previous iteration was a tuning iteration and the current one is not
    // we have reached the end of a tuning phase
    if (_previousIterationWasTuningIteration) {
      ++_numTuningPhasesCompleted;
    }
  }

  _timers.forceUpdateTotal.stop();

  if (potentialEnergy != 0) {
    return std::make_pair(potentialEnergy, virial);
  }
  return std::nullopt;
}

void Simulation::updateVelocities(bool resetForces, RespaIterationType respaIterationType) {
  const double deltaT = _configuration.deltaT.value;

  if (deltaT != 0) {
    _timers.velocityUpdate.start();

    if (respaIterationType == RespaIterationType::OuterStep) {
      TimeDiscretization::calculateVelocities(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                              deltaT, resetForces, /*outerRespaStep*/ true,
                                              _configuration.respaStepSize.value);
    } else {
      TimeDiscretization::calculateVelocities(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                              deltaT, resetForces);
    }

    _timers.velocityUpdate.stop();
  }
}

double Simulation::calculateKineticEnergy() {
  double sum = 0;
  AUTOPAS_OPENMP(parallel shared(_autoPasContainer) reduction(+ : sum))
  for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto mass = _configuration.getParticlePropertiesLibrary()->getMolMass(particle->getTypeId());
    sum += mass * autopas::utils::ArrayMath::dot(particle->getV(), particle->getV());
  }

  return sum * 0.5;
}

void Simulation::updateAngularVelocities(bool resetForces, RespaIterationType respaIterationType,
                                         const size_t outerStepMolID, const size_t innerStepMolID) {
  const double deltaT = _configuration.deltaT.value;

  _timers.angularVelocityUpdate.start();

  if (respaIterationType == RespaIterationType::OuterStep) {
    TimeDiscretization::calculateAngularVelocities(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                                   deltaT, resetForces, /*outerRespaStep*/ true,
                                                   _configuration.respaStepSize.value, outerStepMolID, innerStepMolID);
  } else {
    TimeDiscretization::calculateAngularVelocities(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                                                   deltaT, resetForces, /*outerRespaStep*/ false, 0, outerStepMolID,
                                                   innerStepMolID);
  }

  _timers.angularVelocityUpdate.stop();
}

void Simulation::updateThermostat(bool skipIterationCheck) {
  if (_configuration.useThermostat.value and
      (skipIterationCheck or
       ((not skipIterationCheck) and (_iteration % _configuration.thermostatInterval.value) == 0))) {
    _timers.thermostat.start();
    Thermostat::apply(*_autoPasContainer, *(_configuration.getParticlePropertiesLibrary()),
                      _configuration.targetTemperature.value, _configuration.deltaTemp.value);
    _timers.thermostat.stop();
  }
}

long Simulation::accumulateTime(const long &time) {
  long reducedTime{};
  autopas::AutoPas_MPI_Reduce(&time, &reducedTime, 1, AUTOPAS_MPI_LONG, AUTOPAS_MPI_SUM, 0, AUTOPAS_MPI_COMM_WORLD);

  return reducedTime;
}

std::tuple<bool, std::optional<double>, std::optional<double>> Simulation::calculatePairwiseForces(
    ForceType forceType, bool subtractForces) {
  auto autopasInstanceToUse = _autoPasContainer;
  if (_useSecondAutopasInstanceForRespa) {
    if (forceType == ForceType::CGMolOuter or forceType == ForceType::FPOuter or forceType == ForceType::IBIOuter) {
      autopasInstanceToUse = _autoPasContainerRespa;
    }
  }
  auto wasTuningIterationPotentialEnergyVirial =
      applyWithChosenFunctor<std::tuple<bool, std::optional<double>, std::optional<double>>>(
          [&](auto &&functor) {
            auto wasTuningIteration = autopasInstanceToUse->template computeInteractions(&functor);
            if (functor.getCalculateGlobals()) {
              const auto potentialEnergy = functor.getPotentialEnergy();
              const auto virial = functor.getVirial();
              return std::make_tuple<bool, std::optional<double>, std::optional<double>>(std::move(wasTuningIteration),
                                                                                         potentialEnergy, virial);
            }
            return std::make_tuple<bool, std::optional<double>, std::optional<double>>(std::move(wasTuningIteration),
                                                                                       std::nullopt, std::nullopt);
          },
          forceType, subtractForces);

#if defined(MD_FLEXIBLE_FUNCTOR_COULOMB)
  if (_applyCoulombFunctor and not(forceType == ForceType::CoarseGrain or forceType == ForceType::IBIOuter)) {
    auto wasTuningIterationPotentialEnergyVirialCoulomb =
        applyWithChosenFunctorElectrostatic<std::tuple<bool, std::optional<double>, std::optional<double>>>(
            [&](auto &&functor) {
              auto wasTuningIteration = autopasInstanceToUse->template computeInteractions(&functor);
              if (functor.getCalculateGlobals()) {
                const auto potentialEnergy = functor.getPotentialEnergy();
                const auto virial = functor.getVirial();
                return std::make_tuple<bool, std::optional<double>, std::optional<double>>(
                    std::move(wasTuningIteration), potentialEnergy, virial);
              }
              return std::make_tuple<bool, std::optional<double>, std::optional<double>>(std::move(wasTuningIteration),
                                                                                         std::nullopt, std::nullopt);
            },
            forceType);
    std::get<0>(wasTuningIterationPotentialEnergyVirial) |= std::get<0>(wasTuningIterationPotentialEnergyVirialCoulomb);
    if (std::get<1>(wasTuningIterationPotentialEnergyVirial).has_value()) {
      std::get<1>(wasTuningIterationPotentialEnergyVirial).value() +=
          std::get<1>(wasTuningIterationPotentialEnergyVirialCoulomb).value();
    }
    if (std::get<2>(wasTuningIterationPotentialEnergyVirial).has_value()) {
      std::get<2>(wasTuningIterationPotentialEnergyVirial).value() +=
          std::get<2>(wasTuningIterationPotentialEnergyVirialCoulomb).value();
    }
  }
#endif
  return wasTuningIterationPotentialEnergyVirial;
}

std::tuple<bool, std::optional<double>, std::optional<double>> Simulation::calculateTriwiseForces(ForceType forceType) {
  const auto wasTuningIterationPotentialEnergyVirial =
      applyWithChosenFunctor3B<std::tuple<bool, std::optional<double>, std::optional<double>>>([&](auto &&functor) {
        auto wasTuningIteration = _autoPasContainer->template computeInteractions(&functor);
        if (functor.getCalculateGlobals()) {
          const auto potentialEnergy = functor.getPotentialEnergy();
          const auto virial = functor.getVirial();
          return std::make_tuple<bool, std::optional<double>, std::optional<double>>(std::move(wasTuningIteration),
                                                                                     potentialEnergy, virial);
        }
        return std::make_tuple<bool, std::optional<double>, std::optional<double>>(std::move(wasTuningIteration),
                                                                                   std::nullopt, std::nullopt);
      });
  return wasTuningIterationPotentialEnergyVirial;
}

void Simulation::calculateGlobalForces(const std::array<double, 3> &globalForce) {
  AUTOPAS_OPENMP(parallel shared(_autoPasContainer))
  for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    particle->addF(globalForce);
  }
}

void Simulation::logSimulationState() {
  size_t totalNumberOfParticles{0ul}, ownedParticles{0ul}, haloParticles{0ul},
      totalNumberOfParticlesSecondInstance{0ul}, ownedParticlesSecondInstance{0ul}, haloParticlesSecondInstance{0ul};

  size_t particleCount = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo);
  autopas::AutoPas_MPI_Allreduce(&particleCount, &totalNumberOfParticles, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);

  particleCount = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned);
  autopas::AutoPas_MPI_Allreduce(&particleCount, &ownedParticles, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);

  particleCount = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo);
  autopas::AutoPas_MPI_Allreduce(&particleCount, &haloParticles, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM,
                                 AUTOPAS_MPI_COMM_WORLD);

  if (_useSecondAutopasInstanceForRespa) {
    particleCount = _autoPasContainerRespa->getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo);
    autopas::AutoPas_MPI_Allreduce(&particleCount, &totalNumberOfParticlesSecondInstance, 1, AUTOPAS_MPI_UNSIGNED_LONG,
                                   AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);

    particleCount = _autoPasContainerRespa->getNumberOfParticles(autopas::IteratorBehavior::owned);
    autopas::AutoPas_MPI_Allreduce(&particleCount, &ownedParticlesSecondInstance, 1, AUTOPAS_MPI_UNSIGNED_LONG,
                                   AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);

    particleCount = _autoPasContainerRespa->getNumberOfParticles(autopas::IteratorBehavior::halo);
    autopas::AutoPas_MPI_Allreduce(&particleCount, &haloParticlesSecondInstance, 1, AUTOPAS_MPI_UNSIGNED_LONG,
                                   AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);
  }

  if (_domainDecomposition->getDomainIndex() == 0) {
    std::cout << "\n\n"
              << "Total number of particles at the end of Simulation: " << totalNumberOfParticles << "\n"
              << "Owned: " << ownedParticles << "\n"
              << "Halo : " << haloParticles << "\n";
    if (_useSecondAutopasInstanceForRespa) {
      std::cout << "Total number of particles at the end of Simulation in second Autopas instance: "
                << totalNumberOfParticlesSecondInstance << "\n"
                << "Owned: " << ownedParticlesSecondInstance << "\n"
                << "Halo : " << haloParticlesSecondInstance << "\n";
    }
  }
}

void Simulation::updateSimulationPauseState() {
  // If we are at the beginning of a tuning phase, we need to freeze the simulation
  if (_currentIterationIsTuningIteration and (not _previousIterationWasTuningIteration)) {
    std::cout << "Iteration " << _iteration << ": Freezing simulation for tuning phase. Starting tuning phase...\n";
    _simulationIsPaused = true;
  }

  // If we are at the end of a tuning phase, we need to resume the simulation
  if (_previousIterationWasTuningIteration and (not _currentIterationIsTuningIteration)) {
    std::cout << "Iteration " << _iteration << ": Resuming simulation after tuning phase.\n";

    // reset the forces which accumulated during the tuning phase
    for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      particle->setF(_configuration.globalForce.value);
    }

    // calculate the forces of the latest iteration again
    updateInteractionForces(ForceType::FullParticle);

    _simulationIsPaused = false;
  }
}

void Simulation::logMeasurements() {
  const long positionUpdate = accumulateTime(_timers.positionUpdate.getTotalTime());
  const long quaternionUpdate = accumulateTime(_timers.quaternionUpdate.getTotalTime());
  const long updateContainer = accumulateTime(_timers.updateContainer.getTotalTime());
  const long forceUpdateTotal = accumulateTime(_timers.forceUpdateTotal.getTotalTime());
  const long forceUpdatePairwise = accumulateTime(_timers.forceUpdatePairwise.getTotalTime());
  const long forceUpdateTriwise = accumulateTime(_timers.forceUpdateTriwise.getTotalTime());
  const long forceUpdateTuning = accumulateTime(_timers.forceUpdateTuning.getTotalTime());
  const long forceUpdateNonTuning = accumulateTime(_timers.forceUpdateNonTuning.getTotalTime());
  const long velocityUpdate = accumulateTime(_timers.velocityUpdate.getTotalTime());
  const long angularVelocityUpdate = accumulateTime(_timers.angularVelocityUpdate.getTotalTime());
  const long simulate = accumulateTime(_timers.simulate.getTotalTime());
  const long vtk = accumulateTime(_timers.vtk.getTotalTime());
  const long initialization = accumulateTime(_timers.initialization.getTotalTime());
  const long total = accumulateTime(_timers.total.getTotalTime());
  const long thermostat = accumulateTime(_timers.thermostat.getTotalTime());
  const long haloParticleExchange = accumulateTime(_timers.haloParticleExchange.getTotalTime());
  const long reflectParticlesAtBoundaries = accumulateTime(_timers.reflectParticlesAtBoundaries.getTotalTime());
  const long migratingParticleExchange = accumulateTime(_timers.migratingParticleExchange.getTotalTime());
  const long loadBalancing = accumulateTime(_timers.loadBalancing.getTotalTime());
  const long secondAutopasSync = accumulateTime(_timers.secondAutopasInstanceSync.getTotalTime());

  if (_domainDecomposition->getDomainIndex() == 0) {
    const long wallClockTime = _timers.total.getTotalTime();
    // use the two potentially largest timers to determine the number of chars needed
    const auto maximumNumberOfDigits =
        static_cast<int>(std::max(std::to_string(total).length(), std::to_string(wallClockTime).length()));
    std::cout << "Measurements:\n";
    std::cout << timerToString("Total accumulated                 ", total, maximumNumberOfDigits);
    std::cout << timerToString("  Initialization                  ", initialization, maximumNumberOfDigits, total);
    std::cout << timerToString("  Simulate                        ", simulate, maximumNumberOfDigits, total);
    std::cout << timerToString("    PositionUpdate                ", positionUpdate, maximumNumberOfDigits, simulate);
#if MD_FLEXIBLE_MODE == MULTISITE
    std::cout << timerToString("    QuaternionUpdate              ", quaternionUpdate, maximumNumberOfDigits, simulate);
#endif
    std::cout << timerToString("    UpdateContainer               ", updateContainer, maximumNumberOfDigits, simulate);
    std::cout << timerToString("    Boundaries                    ", haloParticleExchange + migratingParticleExchange,
                               maximumNumberOfDigits, simulate);
    std::cout << timerToString("      HaloParticleExchange        ", haloParticleExchange, maximumNumberOfDigits,
                               haloParticleExchange + reflectParticlesAtBoundaries + migratingParticleExchange);
    std::cout << timerToString("      ReflectParticlesAtBoundaries", reflectParticlesAtBoundaries,
                               maximumNumberOfDigits,
                               haloParticleExchange + reflectParticlesAtBoundaries + migratingParticleExchange);
    std::cout << timerToString("      MigratingParticleExchange   ", migratingParticleExchange, maximumNumberOfDigits,
                               haloParticleExchange + reflectParticlesAtBoundaries + migratingParticleExchange);
    std::cout << timerToString("    ForceUpdateTotal              ", forceUpdateTotal, maximumNumberOfDigits, simulate);
    std::cout << timerToString("      Tuning                      ", forceUpdateTuning, maximumNumberOfDigits,
                               forceUpdateTotal);
    std::cout << timerToString("      ForceUpdateTuning           ", forceUpdateTuning, maximumNumberOfDigits,
                               forceUpdateTotal);
    std::cout << timerToString("      ForceUpdateNonTuning        ", forceUpdateNonTuning, maximumNumberOfDigits,
                               forceUpdateTotal);
    std::cout << timerToString("    SecondAutopasSync             ", secondAutopasSync, maximumNumberOfDigits,
                               simulate);
    std::cout << timerToString("    VelocityUpdate                ", velocityUpdate, maximumNumberOfDigits, simulate);
#if MD_FLEXIBLE_MODE == MULTISITE
    std::cout << timerToString("    AngularVelocityUpdate         ", angularVelocityUpdate, maximumNumberOfDigits,
                               simulate);
#endif
    std::cout << timerToString("    Thermostat                    ", thermostat, maximumNumberOfDigits, simulate);
    std::cout << timerToString("    Vtk                           ", vtk, maximumNumberOfDigits, simulate);
    std::cout << timerToString("    LoadBalancing                 ", loadBalancing, maximumNumberOfDigits, simulate);
    std::cout << timerToString("  One iteration                   ", simulate / static_cast<long>(_iteration),
                               maximumNumberOfDigits, total);

    std::cout << timerToString("Total wall-clock time             ", wallClockTime, maximumNumberOfDigits, total);
    std::cout << "\n";

    std::cout << "Tuning iterations                  : " << _numTuningIterations << " / " << _iteration << " = "
              << (static_cast<double>(_numTuningIterations) / static_cast<double>(_iteration) * 100.) << "%"
              << "\n";

    auto mfups =
        static_cast<double>(_autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned) * _iteration) *
        1e-6 / (static_cast<double>(forceUpdateTotal) * 1e-9);  // 1e-9 for ns to s, 1e-6 for M in MFUPs
    std::cout << "MFUPs/sec                          : " << mfups << "\n";
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    std::cout << "Mean Rebuild Frequency               : " << _autoPasContainer->getMeanRebuildFrequency() << "\n";
#endif
  }
}

bool Simulation::needsMoreIterations() const {
  return _iteration < _configuration.iterations.value or _numTuningPhasesCompleted < _configuration.tuningPhases.value;
}

void Simulation::checkNumParticles(size_t expectedNumParticlesGlobal, size_t numParticlesCurrentlyMigratingLocal,
                                   int lineNumber) {
  if (std::all_of(_configuration.boundaryOption.value.begin(), _configuration.boundaryOption.value.end(),
                  [](const auto &boundary) {
                    return boundary == options::BoundaryTypeOption::periodic or
                           boundary == options::BoundaryTypeOption::reflective;
                  })) {
    const auto numParticlesNowLocal = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned);
    std::array<size_t, 2> sendBuffer{numParticlesNowLocal, numParticlesCurrentlyMigratingLocal};
    std::array<size_t, sendBuffer.size()> receiveBuffer{};
    autopas::AutoPas_MPI_Reduce(sendBuffer.data(), receiveBuffer.data(), receiveBuffer.size(),
                                AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0, AUTOPAS_MPI_COMM_WORLD);
    const auto &[numParticlesNowTotal, numParticlesMigratingTotal] = receiveBuffer;
    if (expectedNumParticlesGlobal != numParticlesNowTotal + numParticlesMigratingTotal) {
      const auto myRank = _domainDecomposition->getDomainIndex();
      std::stringstream ss;
      // clang-format off
      ss << "Rank " << myRank << " Line " << lineNumber
         << ": Particles Lost! All Boundaries are periodic but the number of particles changed:"
         << "Expected        : " << expectedNumParticlesGlobal
         << "Actual          : " << (numParticlesNowTotal + numParticlesMigratingTotal)
         << "  in containers : " << numParticlesNowTotal
         << "  migrating     : " << numParticlesMigratingTotal
         << std::endl;
      // clang-format on
      throw std::runtime_error(ss.str());
    }
  }
}

void Simulation::loadParticles() {
  // Store how many particles are in this config object before removing them.
  const auto numParticlesInConfigLocally = _configuration.particles.size();
  // When loading checkpoints, the config file might contain particles that do not belong to this rank,
  // because rank boundaries don't align with those at the end of the previous simulation due to dynamic load balancing.
  // Idea: Only load what belongs here and send the rest away.
  _autoPasContainer->addParticlesIf(_configuration.particles, [&](auto &p) {
    if (_domainDecomposition->isInsideLocalDomain(p.getR())) {
      // Mark particle in vector as dummy, so we know it has been inserted.
      p.setOwnershipState(autopas::OwnershipState::dummy);
      return true;
    }
    return false;
  });

  // Remove what has been inserted. Everything that remains does not belong into this rank.
  _configuration.particles.erase(std::remove_if(_configuration.particles.begin(), _configuration.particles.end(),
                                                [&](const auto &p) { return p.isDummy(); }),
                                 _configuration.particles.end());

  // Send all remaining particles to all ranks
  // TODO: This is not optimal but since this only happens once upon initialization it is not too bad.
  //       Nevertheless it could be improved by determining which particle has to go to which rank.
  const auto rank = _domainDecomposition->getDomainIndex();
  ParticleCommunicator particleCommunicator(_domainDecomposition->getCommunicator());
  for (int receiverRank = 0; receiverRank < _domainDecomposition->getNumberOfSubdomains(); ++receiverRank) {
    // don't send to ourselves
    if (receiverRank == rank) {
      continue;
    }
    particleCommunicator.sendParticles(_configuration.particles, receiverRank);
  }
  // Erase all locally stored particles. They don't belong to this rank and have been sent away.
  _configuration.flushParticles();

  // Receive particles from all other ranks.
  for (int senderRank = 0; senderRank < _domainDecomposition->getNumberOfSubdomains(); ++senderRank) {
    // don't send to ourselves
    if (senderRank == rank) {
      continue;
    }
    particleCommunicator.receiveParticles(_configuration.particles, senderRank);
  }
  particleCommunicator.waitForSendRequests();

  // Add all new particles that belong in this rank
  _autoPasContainer->addParticlesIf(_configuration.particles,
                                    [&](auto &p) { return _domainDecomposition->isInsideLocalDomain(p.getR()); });
  // cleanup
  _configuration.flushParticles();

  // Output and sanity checks
  const size_t numParticlesLocally = _autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned);
  std::cout << "Number of particles at initialization "
            // align outputs based on the max number of ranks
            << "on rank " << std::setw(std::to_string(_domainDecomposition->getNumberOfSubdomains()).length())
            << std::right << rank << ": " << numParticlesLocally << "\n";

  // Local unnamed struct to pack data for MPI
  struct {
    size_t numParticlesAdded;
    size_t numParticlesInConfig;
  } dataPackage{numParticlesLocally, numParticlesInConfigLocally};
  // Let rank 0 also report the global number of particles
  if (rank == 0) {
    autopas::AutoPas_MPI_Reduce(AUTOPAS_MPI_IN_PLACE, &dataPackage, 2, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0,
                                _domainDecomposition->getCommunicator());
    std::cout << "Number of particles at initialization globally"
              // align ":" with the messages above
              << std::setw(std::to_string(_domainDecomposition->getNumberOfSubdomains()).length()) << ""
              << ": " << dataPackage.numParticlesAdded << "\n";
    // Sanity check that on a global scope all particles have been loaded
    if (dataPackage.numParticlesAdded != dataPackage.numParticlesInConfig) {
      throw std::runtime_error(
          "Simulation::loadParticles(): "
          "Not all particles from the configuration file could be added to AutoPas!\n"
          "Configuration : " +
          std::to_string(dataPackage.numParticlesInConfig) +
          "\n"
          "Added globally: " +
          std::to_string(dataPackage.numParticlesAdded));
    }
  } else {
    // In-place reduce needs different calls on root vs rest...
    autopas::AutoPas_MPI_Reduce(&dataPackage, nullptr, 2, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0,
                                _domainDecomposition->getCommunicator());
  }
}

template <class ReturnType, class FunctionType>
ReturnType Simulation::applyWithChosenFunctor(FunctionType f, ForceType forceType, bool subtractForces) {
  const bool useCGFunctor = forceType == ForceType::CoarseGrain;
  const double cutoff = _configuration.cutoff.value;
  auto &particlePropertiesLibrary = *_configuration.getParticlePropertiesLibrary();
  if (useCGFunctor) {
    if (_configuration.useApproxForceRespa.value) {
      auto func = LJFunctorTypeAutovecApproxMultisite{cutoff, particlePropertiesLibrary, subtractForces ? -1 : 1};
      return f(func);
    } else {
      auto func = LuTFunctorType{cutoff, particlePropertiesLibrary, subtractForces ? -1 : 1};
      func.setLuT(_lut);
      return f(func);
    }
  }

  if (forceType == ForceType::IBIOuter) {
    auto func = LuTFunctorType{cutoff * _configuration.cutoffFactorRespa.value, particlePropertiesLibrary,
                               subtractForces ? -1 : 1};
    func.setInnerCutoff(cutoff);
    func.setLuT(_lut);
    return f(func);
  }

  switch (_configuration.functorOption.value) {
    case MDFlexConfig::FunctorOption::lj12_6: {
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
      if (forceType == ForceType::FPOuter or forceType == ForceType::CGMolOuter) {
        auto func = LJFunctorTypeAutovec{cutoff * _configuration.cutoffFactorRespa.value, particlePropertiesLibrary};
        func.setInnerCutoff(cutoff);
        return f(func);
      } else {
        return f(LJFunctorTypeAutovec{cutoff, particlePropertiesLibrary});
      }

#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor AutoVec. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AUTOVEC=ON`.");
#endif
    }
    case MDFlexConfig::FunctorOption::lj12_6_AVX: {
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
      return f(LJFunctorTypeAVX{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor AVX. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AVX=ON`.");
#endif
    }
    case MDFlexConfig::FunctorOption::lj12_6_SVE: {
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
      return f(LJFunctorTypeSVE{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor SVE. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_SVE=ON`.");
#endif
    }
    default: {
      throw std::runtime_error("Unknown pairwise functor choice" +
                               std::to_string(static_cast<int>(_configuration.functorOption.value)));
    }
  }
}

template <class ReturnType, class FunctionType>
ReturnType Simulation::applyWithChosenFunctorElectrostatic(FunctionType f, ForceType forceType) {
  const double cutoff = _configuration.cutoff.value;
  auto &particlePropertiesLibrary = *_configuration.getParticlePropertiesLibrary();

#if defined(MD_FLEXIBLE_FUNCTOR_COULOMB)
  if (forceType == ForceType::FPOuter or forceType == ForceType::CGMolOuter) {
    auto func = CoulombFunctorTypeAutovec{cutoff * _configuration.cutoffFactorRespa.value, particlePropertiesLibrary};
    func.setInnerCutoff(cutoff);
    return f(func);
  } else {
    return f(CoulombFunctorTypeAutovec{cutoff, particlePropertiesLibrary});
  }
#else
  throw std::runtime_error(
      "MD-Flexible was not compiled with support for Coulomb interactions. Activate it via `cmake "
      "-DMD_FLEXIBLE_FUNCTOR_COULOMB=ON`.");
#endif
}

template <class ReturnType, class FunctionType>
ReturnType Simulation::applyWithChosenFunctor3B(FunctionType f, bool useCGFunctor) {
  const double cutoff = _configuration.cutoff.value;
  auto &particlePropertiesLibrary = *_configuration.getParticlePropertiesLibrary();
  switch (_configuration.functorOption3B.value) {
    case MDFlexConfig::FunctorOption3B::at: {
#if defined(MD_FLEXIBLE_FUNCTOR_AT_AUTOVEC)
      return f(ATFunctor{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for AxilrodTeller Functor. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AT_AUTOVEC=ON`.");
#endif
    }
    default: {
      throw std::runtime_error("Unknown triwise functor choice" +
                               std::to_string(static_cast<int>(_configuration.functorOption3B.value)));
    }
  }
}
