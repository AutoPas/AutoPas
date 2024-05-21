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
#include <vector>

// The LogicHandler includes dependencies to wide parts of AutoPas, making it expensive to compile and thus is moved
// here from AutoPasDecl.h.
#include "autopas/AutoPasDecl.h"
#include "autopas/InstanceCounter.h"
#include "autopas/LogicHandler.h"
#include "autopas/Version.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyFactory.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyLogger.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/CompileInfo.h"
#include "autopas/utils/NumberInterval.h"
#include "autopas/utils/NumberSetFinite.h"
#include "autopas/utils/WrapOpenMP.h"

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
  _autoTuners = std::move(other._autoTuners);
  _logicHandler = std::move(other._logicHandler);
  return *this;
}

template <class Particle>
void AutoPas<Particle>::init() {
  AutoPasLog(INFO, "AutoPas Version: {}", AutoPas_VERSION);
  AutoPasLog(INFO, "Compiled with  : {}", utils::CompileInfo::getCompilerInfo());

  if (_tuningStrategyFactoryInfo.autopasMpiCommunicator == AUTOPAS_MPI_COMM_NULL) {
    AutoPas_MPI_Comm_dup(AUTOPAS_MPI_COMM_WORLD, &_tuningStrategyFactoryInfo.autopasMpiCommunicator);
  } else {
    _externalMPICommunicator = true;
  }
  if (std::find(_tuningStrategyOptions.begin(), _tuningStrategyOptions.end(),
                TuningStrategyOption::mpiDivideAndConquer) != _tuningStrategyOptions.end()) {
    _tuningStrategyFactoryInfo.mpiDivideAndConquer = true;
  }

  _logicHandlerInfo.sortingThreshold = _sortingThreshold;

  // If an interval was given for the cell size factor, change it to the relevant values.
  // Don't modify _allowedCellSizeFactors to preserve the initial (type) information.
  const auto cellSizeFactors = [&]() -> NumberSetFinite<double> {
    if (const auto *csfIntervalPtr = dynamic_cast<NumberInterval<double> *>(_allowedCellSizeFactors.get())) {
      const auto interactionLength =
          _logicHandlerInfo.cutoff * _logicHandlerInfo.verletSkinPerTimestep * _verletRebuildFrequency;
      const auto boxLengthX = _logicHandlerInfo.boxMax[0] - _logicHandlerInfo.boxMin[0];
      return {SearchSpaceGenerators::calculateRelevantCsfs(*csfIntervalPtr, interactionLength, boxLengthX)};
    } else {
      // in this case _allowedCellSizeFactors is a finite set
      return {_allowedCellSizeFactors->getAll()};
    }
  }();

  // Create autotuners for each interaction type
  for (const auto &interactionType : _allowedInteractionTypeOptions) {
    const auto searchSpace = SearchSpaceGenerators::cartesianProduct(
        _allowedContainers, _allowedTraversals[interactionType], _allowedLoadEstimators,
        _allowedDataLayouts[interactionType], _allowedNewton3Options[interactionType], &cellSizeFactors,
        interactionType);

    AutoTuner::TuningStrategiesListType tuningStrategies;
    tuningStrategies.reserve(_tuningStrategyOptions.size());
    for (const auto &strategy : _tuningStrategyOptions) {
      tuningStrategies.emplace_back(TuningStrategyFactory::generateTuningStrategy(
          searchSpace, strategy, _tuningStrategyFactoryInfo, _outputSuffix));
    }
    if (_useTuningStrategyLoggerProxy) {
      tuningStrategies.emplace_back(std::make_unique<TuningStrategyLogger>(_outputSuffix));
    }
    _autoTuners.emplace(interactionType,
                        std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, _autoTunerInfo,
                                                             _verletRebuildFrequency, _outputSuffix));
  }

  // Create logic handler
  _logicHandler = std::make_unique<std::remove_reference_t<decltype(*_logicHandler)>>(
      _autoTuners, _logicHandlerInfo, _verletRebuildFrequency, _outputSuffix);
}

template <class Particle>
template <class Functor>
bool AutoPas<Particle>::computeInteractions(Functor *f) {
  static_assert(
      not std::is_same_v<Functor, autopas::Functor<Particle, Functor>>,
      "The static type of Functor in computeInteractions is not allowed to be autopas::Functor. Please use the "
      "derived type instead, e.g. by using a dynamic_cast.");
  if (f->getCutoff() > this->getCutoff()) {
    utils::ExceptionHandler::exception("Functor cutoff ({}) must not be larger than container cutoff ({})",
                                       f->getCutoff(), this->getCutoff());
  }

  if constexpr (utils::isPairwiseFunctor<Functor>()) {
    return _logicHandler->template computeInteractionsPipeline<Functor>(f, InteractionTypeOption::pairwise);
  } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
    return _logicHandler->template computeInteractionsPipeline<Functor>(f, InteractionTypeOption::triwise);
  } else {
    utils::ExceptionHandler::exception(
        "Functor is not valid. Only 2-body and 3-body functors are supported. Please use a functor derived from "
        "PairwiseFunctor or TriwiseFunctor.");
  }
  return false;
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
  AUTOPAS_OPENMP(parallel for schedule(static, std::max(1ul, collectionSize / omp_get_max_threads())))
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
  AUTOPAS_OPENMP(parallel for reduction(+ : numTrue))
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
  if (_allowedCellSizeFactors->isInterval()) {
    AutoPasLog(WARN,
               "The allowed Cell Size Factors are a continuous interval but internally only those values that "
               "yield unique numbers of cells are used. Resizing does not cause these values to be recalculated so "
               "the same configurations might now yield different and non-unique numbers of cells!");
  }
  _logicHandlerInfo.boxMin = boxMin;
  _logicHandlerInfo.boxMax = boxMax;
  return _logicHandler->resizeBox(boxMin, boxMax);
}

template <class Particle>
void AutoPas<Particle>::forceRetune() {
  for (auto &[interaction, tuner] : _autoTuners) {
    tuner->forceRetune();
  }
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
  AUTOPAS_OPENMP(parallel for reduction(+ : numTrue))
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
  bool isTrivial = true;
  for (auto &[interaction, tuner] : _autoTuners) {
    isTrivial = isTrivial and tuner->searchSpaceIsTrivial();
  }
  return isTrivial;
}

}  // namespace autopas
