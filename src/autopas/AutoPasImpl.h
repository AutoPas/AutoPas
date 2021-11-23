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
#include "autopas/Version.h"

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
  AutoPasLog(info, "AutoPas Version: {}", AutoPas_VERSION);
  if (_numSamples % _verletRebuildFrequency != 0) {
    AutoPasLog(warn,
               "Number of samples ({}) is not a multiple of the rebuild frequency ({}). This can lead to problems "
               "when multiple AutoPas instances interact (e.g. via MPI).",
               _numSamples, _verletRebuildFrequency);
  }

  if (_autopasMPICommunicator == AUTOPAS_MPI_COMM_NULL) {
    AutoPas_MPI_Comm_dup(AUTOPAS_MPI_COMM_WORLD, &_autopasMPICommunicator);
  } else {
    _externalMPICommunicator = true;
  }
  _autoTuner = std::make_unique<autopas::AutoTuner<Particle>>(
      _boxMin, _boxMax, _cutoff, _verletSkin, _verletClusterSize,
      std::move(TuningStrategyFactory::generateTuningStrategy(
          _tuningStrategyOption, _allowedContainers, *_allowedCellSizeFactors, _allowedTraversals,
          _allowedLoadEstimators, _allowedDataLayouts, _allowedNewton3Options, _maxEvidence, _relativeOptimumRange,
          _maxTuningPhasesWithoutTest, _relativeBlacklistRange, _evidenceFirstPrediction, _acquisitionFunctionOption,
          _extrapolationMethodOption, _outputSuffix, _mpiStrategyOption, _autopasMPICommunicator)),
      _selectorStrategy, _tuningMetricOption, _tuningInterval, _numSamples, _outputSuffix);
  _logicHandler =
      std::make_unique<std::remove_reference_t<decltype(*_logicHandler)>>(*(_autoTuner.get()), _verletRebuildFrequency);
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
  return _logicHandler->iteratePairwise(f);
}

template <class Particle>
size_t AutoPas<Particle>::getNumberOfParticles(IteratorBehavior behavior) const {
  size_t numParticles{0};
  if (behavior & IteratorBehavior::owned) {
    numParticles += _logicHandler->getNumParticlesOwned();
  }
  if (behavior & IteratorBehavior::halo) {
    numParticles += _logicHandler->getNumParticlesHalo();
  }
  // non fatal sanity check whether the behavior contained anything else
  if (behavior & ~(IteratorBehavior::ownedOrHalo)) {
    utils::ExceptionHandler::exception(
        "AutoPas::getNumberOfParticles() does not support iterator behaviors other than owned or halo.");
  }

  return numParticles;
}

template <class Particle>
void AutoPas<Particle>::addParticle(const Particle &p) {
  _logicHandler->addParticle(p);
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
void AutoPas<Particle>::deleteAllParticles() {
  _logicHandler->deleteAllParticles();
}

template <class Particle>
void AutoPas<Particle>::deleteParticle(ParticleIteratorWrapper<Particle, true> &iter) {
  _logicHandler->deleteParticle(iter);
}

template <class Particle>
typename AutoPas<Particle>::iterator_t AutoPas<Particle>::begin(IteratorBehavior behavior) {
  return _logicHandler->begin(behavior);
}

template <class Particle>
typename AutoPas<Particle>::const_iterator_t AutoPas<Particle>::begin(IteratorBehavior behavior) const {
  return std::as_const(*_logicHandler).begin(behavior);
}

template <class Particle>
typename AutoPas<Particle>::iterator_t AutoPas<Particle>::getRegionIterator(std::array<double, 3> lowerCorner,
                                                                            std::array<double, 3> higherCorner,
                                                                            IteratorBehavior behavior) {
  return _logicHandler->getRegionIterator(lowerCorner, higherCorner, behavior);
}

template <class Particle>
typename AutoPas<Particle>::const_iterator_t AutoPas<Particle>::getRegionIterator(std::array<double, 3> lowerCorner,
                                                                                  std::array<double, 3> higherCorner,
                                                                                  IteratorBehavior behavior) const {
  return std::as_const(*_logicHandler).getRegionIterator(lowerCorner, higherCorner, behavior);
}

template <class Particle>
unsigned long AutoPas<Particle>::getContainerType() const {
  return _autoTuner->getContainer()->getContainerType();
}

template <class Particle>
std::array<double, 3> AutoPas<Particle>::getBoxMin() const {
  return _autoTuner->getContainer()->getBoxMin();
}

template <class Particle>
std::array<double, 3> AutoPas<Particle>::getBoxMax() const {
  return _autoTuner->getContainer()->getBoxMax();
}

template <class Particle>
std::shared_ptr<autopas::ParticleContainerInterface<Particle>> AutoPas<Particle>::getContainer() {
  return _autoTuner->getContainer();
}

template <class Particle>
const std::shared_ptr<autopas::ParticleContainerInterface<Particle>> AutoPas<Particle>::getContainer() const {
  return _autoTuner->getContainer();
}

}  // namespace autopas
