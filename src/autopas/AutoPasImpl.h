/**
 * @file AutoPasImpl.h
 * Contains all the missing implementations from AutoPasDecl.h that cannot be implemented there because AutoTuner and
 * LogicHandler are only forward declared.
 */

#pragma once

#include <array>
#include <memory>
#include <ostream>
#include <type_traits>

#include "autopas/AutoPasDecl.h"
#include "autopas/InstanceCounter.h"
#include "autopas/LogicHandlerInfo.h"
#include "autopas/Version.h"
#include "autopas/utils/CompileInfo.h"

// These next three includes have dependencies to all of AutoPas and thus are moved here from AutoPasDecl.h.
#include "autopas/LogicHandler.h"
#include "autopas/selectors/AutoTuner.h"
#include "autopas/selectors/tuningStrategy/TuningStrategyFactory.h"

namespace autopas {

template <class Particle>
AutoPas<Particle>::AutoPas(std::ostream &logOutputStream) {
  // count the number of autopas instances. This is needed to ensure that the autopas
  // logger is not unregistered while other instances are still using it.
  InstanceCounter::count++;
  // remove potentially existing logger
  autopas::Logger::unregister();
  // initialize the Logger
  autopas::Logger::create(logOutputStream);
  // The logger is normally only flushed on successful program termination.
  // This line ensures flushing when log messages of level warning or more severe are created.
  autopas::Logger::get()->flush_on(spdlog::level::warn);
}

template <class Particle>
AutoPas<Particle>::~AutoPas() {
  InstanceCounter::count--;
  if (InstanceCounter::count == 0) {
    // remove the Logger from the registry. Do this only if we have no other autopas instances running.
    autopas::Logger::unregister();
  }
}

template <class Particle>
AutoPas<Particle> &AutoPas<Particle>::operator=(AutoPas &&other) noexcept {
  _autoTuner = std::move(other._autoTuner);
  _logicHandler = std::move(other._logicHandler);
  return *this;
}

template <class Particle>
void AutoPas<Particle>::init() {
  AutoPasLog(INFO, "AutoPas Version: {}", AutoPas_VERSION);
  AutoPasLog(INFO, "Compiled with  : {}", utils::CompileInfo::getCompilerInfo());
  if (_numSamples % _verletRebuildFrequency != 0) {
    AutoPasLog(WARN,
               "Number of samples ({}) is not a multiple of the rebuild frequency ({}). This can lead to problems "
               "when multiple AutoPas instances interact (e.g. via MPI).",
               _numSamples, _verletRebuildFrequency);
  }

  if (_autopasMPICommunicator == AUTOPAS_MPI_COMM_NULL) {
    AutoPas_MPI_Comm_dup(AUTOPAS_MPI_COMM_WORLD, &_autopasMPICommunicator);
  } else {
    _externalMPICommunicator = true;
  }
  auto tuningStrategy = TuningStrategyFactory::generateTuningStrategy(
      _tuningStrategyOption, _allowedContainers, *_allowedCellSizeFactors, _allowedTraversals, _allowedLoadEstimators,
      _allowedDataLayouts, _allowedNewton3Options, _maxEvidence, _relativeOptimumRange, _maxTuningPhasesWithoutTest,
      _relativeBlacklistRange, _evidenceFirstPrediction, _acquisitionFunctionOption, _extrapolationMethodOption,
      _ruleFileName, _outputSuffix, _mpiStrategyOption, _autopasMPICommunicator);
  _autoTuner = std::make_unique<autopas::AutoTuner<Particle>>(
      std::move(tuningStrategy), _mpiTuningMaxDifferenceForBucket, _mpiTuningWeightForMaxDensity, _selectorStrategy,
      _tuningMetricOption, _tuningInterval, _numSamples, _verletRebuildFrequency, _outputSuffix, _useTuningLogger);
  LogicHandlerInfo logicHandlerInfo{
      _boxMin, _boxMax, _cutoff, _verletSkinPerTimestep, _verletRebuildFrequency, _verletClusterSize, _outputSuffix};
  _logicHandler = std::make_unique<std::remove_reference_t<decltype(*_logicHandler)>>(*_autoTuner, logicHandlerInfo);
}

template <class Particle>
template <class Functor>
bool AutoPas<Particle>::iteratePairwise(Functor *f) {
  static_assert(not std::is_same<Functor, autopas::Functor<Particle, Functor>>::value,
                "The static type of Functor in iteratePairwise is not allowed to be autopas::Functor. Please use the "
                "derived type instead, e.g. by using a dynamic_cast.");
  if (f->getCutoff() > this->getCutoff()) {
    utils::ExceptionHandler::exception("Functor cutoff ({}) must not be larger than container cutoff ({})",
                                       f->getCutoff(), this->getCutoff());
  }
  return _logicHandler->iteratePairwisePipeline(f);
}

template <class Particle>
size_t AutoPas<Particle>::getNumberOfParticles(IteratorBehavior behavior) const {
  size_t numParticles{0};
  if (behavior & IteratorBehavior::owned) {
    numParticles += _logicHandler->getNumberOfParticlesOwned();
  }
  if (behavior & IteratorBehavior::halo) {
    numParticles += _logicHandler->getNumberOfParticlesHalo();
  }
  // non fatal sanity check whether the behavior contained anything else
  if (behavior & ~(IteratorBehavior::ownedOrHalo)) {
    utils::ExceptionHandler::exception(
        "AutoPas::getNumberOfParticles() does not support iterator behaviors other than owned or halo.");
  }

  return numParticles;
}

template <class Particle>
void AutoPas<Particle>::reserve(size_t numParticles) {
  _logicHandler->reserve(numParticles);
}

template <class Particle>
void AutoPas<Particle>::reserve(size_t numParticles, size_t numHaloParticles) {
  _logicHandler->reserve(numParticles, numHaloParticles);
}

template <class Particle>
template <class F>
void AutoPas<Particle>::addParticlesAux(size_t numParticlesToAdd, size_t numHalosToAdd, size_t collectionSize,
                                        F loopBody) {
  reserve(getNumberOfParticles(IteratorBehavior::owned) + numParticlesToAdd,
          getNumberOfParticles(IteratorBehavior::halo) + numHalosToAdd);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for schedule(static, std::max(1ul, collectionSize / omp_get_max_threads()))
#endif
  for (auto i = 0; i < collectionSize; ++i) {
    loopBody(i);
  }
}

template <class Particle>
void AutoPas<Particle>::addParticle(const Particle &p) {
  _logicHandler->addParticle(p);
}

template <class Particle>
template <class Collection>
void AutoPas<Particle>::addParticles(Collection &&particles) {
  addParticlesAux(particles.size(), 0, particles.size(), [&](auto i) { addParticle(particles[i]); });
}

template <class Particle>
template <class Collection, class F>
void AutoPas<Particle>::addParticlesIf(Collection &&particles, F predicate) {
  std::vector<char> predicateMask(particles.size());
  int numTrue = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(+ : numTrue)
#endif
  for (auto i = 0; i < particles.size(); ++i) {
    if (predicate(particles[i])) {
      predicateMask[i] = static_cast<char>(true);
      ++numTrue;
    } else {
      predicateMask[i] = static_cast<char>(false);
    }
  }

  addParticlesAux(numTrue, 0, particles.size(), [&](auto i) {
    if (predicateMask[i]) {
      addParticle(particles[i]);
    }
  });
}

template <class Particle>
std::vector<Particle> AutoPas<Particle>::updateContainer() {
  return _logicHandler->updateContainer();
}

template <class Particle>
std::vector<Particle> AutoPas<Particle>::resizeBox(const std::array<double, 3> &boxMin,
                                                   const std::array<double, 3> &boxMax) {
  _boxMin = boxMin;
  _boxMax = boxMax;
  return _logicHandler->resizeBox(boxMin, boxMax);
}

template <class Particle>
void AutoPas<Particle>::forceRetune() {
  _autoTuner->forceRetune();
}

template <class Particle>
void AutoPas<Particle>::addHaloParticle(const Particle &haloParticle) {
  _logicHandler->addHaloParticle(haloParticle);
}

template <class Particle>
template <class Collection>
void AutoPas<Particle>::addHaloParticles(Collection &&particles) {
  addParticlesAux(0, particles.size(), particles.size(), [&](auto i) { addHaloParticle(particles[i]); });
}

template <class Particle>
template <class Collection, class F>
void AutoPas<Particle>::addHaloParticlesIf(Collection &&particles, F predicate) {
  std::vector<char> predicateMask(particles.size());
  int numTrue = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(+ : numTrue)
#endif
  for (auto i = 0; i < particles.size(); ++i) {
    if (predicate(particles[i])) {
      predicateMask[i] = static_cast<char>(true);
      ++numTrue;
    } else {
      predicateMask[i] = static_cast<char>(false);
    }
  }

  addParticlesAux(0, numTrue, particles.size(), [&](auto i) {
    if (predicateMask[i]) {
      addHaloParticle(particles[i]);
    }
  });
}

template <class Particle>
void AutoPas<Particle>::deleteAllParticles() {
  _logicHandler->deleteAllParticles();
}

template <class Particle>
void AutoPas<Particle>::deleteParticle(IteratorT &iter) {
  _logicHandler->decreaseParticleCounter(*iter);
  internal::deleteParticle(iter);
}

template <class Particle>
void AutoPas<Particle>::deleteParticle(RegionIteratorT &iter) {
  _logicHandler->decreaseParticleCounter(*iter);
  internal::deleteParticle(iter);
}

template <class Particle>
bool AutoPas<Particle>::deleteParticle(Particle &particle) {
  _logicHandler->decreaseParticleCounter(particle);
  // if the particle was not found in the logic handler's buffers it must be in the container
  auto [particleDeleted, refStillValid] = _logicHandler->deleteParticleFromBuffers(particle);
  if (not particleDeleted) {
    refStillValid = _logicHandler->getContainer().deleteParticle(particle);
  }
  return refStillValid;
}

template <class Particle>
typename AutoPas<Particle>::IteratorT AutoPas<Particle>::begin(IteratorBehavior behavior) {
  return _logicHandler->begin(behavior);
}

template <class Particle>
typename AutoPas<Particle>::ConstIteratorT AutoPas<Particle>::begin(IteratorBehavior behavior) const {
  return std::as_const(*_logicHandler).begin(behavior);
}

template <class Particle>
typename AutoPas<Particle>::RegionIteratorT AutoPas<Particle>::getRegionIterator(
    const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
  return _logicHandler->getRegionIterator(lowerCorner, higherCorner, behavior);
}

template <class Particle>
typename AutoPas<Particle>::RegionConstIteratorT AutoPas<Particle>::getRegionIterator(
    const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
    IteratorBehavior behavior) const {
  return std::as_const(*_logicHandler).getRegionIterator(lowerCorner, higherCorner, behavior);
}

template <class Particle>
unsigned long AutoPas<Particle>::getContainerType() const {
  return _logicHandler->getContainer().getContainerType();
}

template <class Particle>
const std::array<double, 3> &AutoPas<Particle>::getBoxMin() const {
  return _logicHandler->getContainer().getBoxMin();
}

template <class Particle>
const std::array<double, 3> &AutoPas<Particle>::getBoxMax() const {
  return _logicHandler->getContainer().getBoxMax();
}

template <class Particle>
autopas::ParticleContainerInterface<Particle> &AutoPas<Particle>::getContainer() {
  return _logicHandler->getContainer();
}

template <class Particle>
const autopas::ParticleContainerInterface<Particle> &AutoPas<Particle>::getContainer() const {
  return _logicHandler->getContainer();
}

template <class Particle>
bool AutoPas<Particle>::searchSpaceIsTrivial() {
  return _autoTuner->searchSpaceIsTrivial();
}

}  // namespace autopas
