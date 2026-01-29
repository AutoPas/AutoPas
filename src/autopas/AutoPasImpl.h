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

template <class Particle_T>
AutoPas<Particle_T>::AutoPas(std::ostream &logOutputStream) {
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

template <class Particle_T>
AutoPas<Particle_T>::~AutoPas() {
  InstanceCounter::count--;
  if (InstanceCounter::count == 0) {
    // remove the Logger from the registry. Do this only if we have no other autopas instances running.
    autopas::Logger::unregister();
  }
}

template <class Particle_T>
AutoPas<Particle_T> &AutoPas<Particle_T>::operator=(AutoPas &&other) noexcept {
  _autoTuners = std::move(other._autoTuners);
  _logicHandler = std::move(other._logicHandler);
  return *this;
}

template <class Particle_T>
void AutoPas<Particle_T>::init() {
  int myRank{};
  AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  if (myRank == 0) {
    AutoPasLog(INFO, "AutoPas Version: {}", AutoPas_VERSION);
    AutoPasLog(INFO, "Compiled with  : {}", utils::CompileInfo::getCompilerInfo());
  }

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
  _logicHandlerInfo.useMortonIndex = _useMortonIndex;
  _logicHandlerInfo.useLiveId = _useLiveId;
  _logicHandlerInfo.preloadLJMixingPtr = _preloadLJMixingPtr;
  _logicHandlerInfo.reserveVLSizes = _reserveVLSizes;

  // If an interval was given for the cell size factor, change it to the relevant values.
  // Don't modify _allowedCellSizeFactors to preserve the initial (type) information.
  const auto cellSizeFactors = [&]() -> NumberSetFinite<double> {
    if (const auto *csfIntervalPtr = dynamic_cast<NumberInterval<double> *>(_allowedCellSizeFactors.get())) {
      const auto interactionLength = _logicHandlerInfo.cutoff * _logicHandlerInfo.verletSkin;
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
          searchSpace, strategy, _tuningStrategyFactoryInfo, interactionType, _outputSuffix));
    }
    if (_useTuningStrategyLoggerProxy) {
      tuningStrategies.emplace_back(std::make_unique<TuningStrategyLogger>(_outputSuffix));
    }
    auto tunerOutputSuffix = _outputSuffix + "_" + interactionType.to_string();
    _autoTuners.emplace(interactionType,
                        std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, _autoTunerInfo,
                                                             _verletRebuildFrequency, tunerOutputSuffix));
  }

  // Create logic handler
  _logicHandler = std::make_unique<std::remove_reference_t<decltype(*_logicHandler)>>(
      _autoTuners, _logicHandlerInfo, _verletRebuildFrequency, _outputSuffix);
}

template <class Particle_T>
template <class Functor>
bool AutoPas<Particle_T>::computeInteractions(Functor *f) {
  static_assert(
      not std::is_same_v<Functor, autopas::Functor<Particle_T, Functor>>,
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
        "Functor is not valid. Only pairwise and triwise functors are supported. Please use a functor derived from "
        "PairwiseFunctor or TriwiseFunctor.");
  }
  return false;
}

template <class Particle_T>
size_t AutoPas<Particle_T>::getNumberOfParticles(IteratorBehavior behavior) const {
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

template <class Particle_T>
void AutoPas<Particle_T>::reserve(size_t numParticles) {
  _logicHandler->reserve(numParticles);
}

template <class Particle_T>
void AutoPas<Particle_T>::reserve(size_t numParticles, size_t numHaloParticles) {
  _logicHandler->reserve(numParticles, numHaloParticles);
}

template <class Particle_T>
template <class F>
void AutoPas<Particle_T>::addParticlesAux(size_t numParticlesToAdd, size_t numHalosToAdd, size_t collectionSize,
                                          F loopBody) {
  reserve(getNumberOfParticles(IteratorBehavior::owned) + numParticlesToAdd,
          getNumberOfParticles(IteratorBehavior::halo) + numHalosToAdd);
  AUTOPAS_OPENMP(parallel for schedule(static, std::max(1ul, collectionSize / omp_get_max_threads())))
  for (auto i = 0; i < collectionSize; ++i) {
    loopBody(i);
  }
}

template <class Particle_T>
void AutoPas<Particle_T>::addParticle(const Particle_T &p) {
  _logicHandler->addParticle(p);
}

template <class Particle_T>
template <class Collection>
void AutoPas<Particle_T>::addParticles(Collection &&particles) {
  addParticlesAux(particles.size(), 0, particles.size(), [&](auto i) { addParticle(particles[i]); });
}

template <class Particle_T>
template <class Collection, class F>
void AutoPas<Particle_T>::addParticlesIf(Collection &&particles, F predicate) {
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

template <class Particle_T>
std::vector<Particle_T> AutoPas<Particle_T>::updateContainer() {
  return _logicHandler->updateContainer();
}

template <class Particle_T>
std::vector<Particle_T> AutoPas<Particle_T>::resizeBox(const std::array<double, 3> &boxMin,
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

template <class Particle_T>
void AutoPas<Particle_T>::forceRetune() {
  for (auto &[interaction, tuner] : _autoTuners) {
    tuner->forceRetune();
  }
}

template <class Particle_T>
void AutoPas<Particle_T>::addHaloParticle(const Particle_T &haloParticle) {
  _logicHandler->addHaloParticle(haloParticle);
}

template <class Particle_T>
template <class Collection>
void AutoPas<Particle_T>::addHaloParticles(Collection &&particles) {
  addParticlesAux(0, particles.size(), particles.size(), [&](auto i) { addHaloParticle(particles[i]); });
}

template <class Particle_T>
template <class Collection, class F>
void AutoPas<Particle_T>::addHaloParticlesIf(Collection &&particles, F predicate) {
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

template <class Particle_T>
void AutoPas<Particle_T>::deleteAllParticles() {
  _logicHandler->deleteAllParticles();
}

template <class Particle_T>
void AutoPas<Particle_T>::deleteParticle(IteratorT &iter) {
  _logicHandler->decreaseParticleCounter(*iter);
  internal::deleteParticle(iter);
}

template <class Particle_T>
void AutoPas<Particle_T>::deleteParticle(RegionIteratorT &iter) {
  _logicHandler->decreaseParticleCounter(*iter);
  internal::deleteParticle(iter);
}

template <class Particle_T>
bool AutoPas<Particle_T>::deleteParticle(Particle_T &particle) {
  _logicHandler->decreaseParticleCounter(particle);
  // if the particle was not found in the logic handler's buffers it must be in the container
  auto [particleDeleted, refStillValid] = _logicHandler->deleteParticleFromBuffers(particle);
  if (not particleDeleted) {
    refStillValid = _logicHandler->getContainer().deleteParticle(particle);
  }
  return refStillValid;
}

template <class Particle_T>
typename AutoPas<Particle_T>::IteratorT AutoPas<Particle_T>::begin(IteratorBehavior behavior) {
  return _logicHandler->begin(behavior);
}

template <class Particle_T>
typename AutoPas<Particle_T>::ConstIteratorT AutoPas<Particle_T>::begin(IteratorBehavior behavior) const {
  return std::as_const(*_logicHandler).begin(behavior);
}

template <class Particle_T>
typename AutoPas<Particle_T>::RegionIteratorT AutoPas<Particle_T>::getRegionIterator(
    const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
  return _logicHandler->getRegionIterator(lowerCorner, higherCorner, behavior);
}

template <class Particle_T>
typename AutoPas<Particle_T>::RegionConstIteratorT AutoPas<Particle_T>::getRegionIterator(
    const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
    IteratorBehavior behavior) const {
  return std::as_const(*_logicHandler).getRegionIterator(lowerCorner, higherCorner, behavior);
}

template <class Particle_T>
unsigned long AutoPas<Particle_T>::getContainerType() const {
  return _logicHandler->getContainer().getContainerType();
}

template <class Particle_T>
const std::array<double, 3> &AutoPas<Particle_T>::getBoxMin() const {
  return _logicHandler->getContainer().getBoxMin();
}

template <class Particle_T>
const std::array<double, 3> &AutoPas<Particle_T>::getBoxMax() const {
  return _logicHandler->getContainer().getBoxMax();
}

template <class Particle_T>
autopas::ParticleContainerInterface<Particle_T> &AutoPas<Particle_T>::getContainer() {
  return _logicHandler->getContainer();
}

template <class Particle_T>
const autopas::ParticleContainerInterface<Particle_T> &AutoPas<Particle_T>::getContainer() const {
  return _logicHandler->getContainer();
}

template <class Particle_T>
bool AutoPas<Particle_T>::searchSpaceIsTrivial() {
  bool isTrivial = true;
  for (auto &[interaction, tuner] : _autoTuners) {
    isTrivial = isTrivial and tuner->searchSpaceIsTrivial();
  }
  return isTrivial;
}

}  // namespace autopas
