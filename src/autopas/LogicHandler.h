/**
 * @file LogicHandler.h
 * @author seckler
 * @date 31.05.19
 */

#pragma once
#include <limits>

#include "autopas/iterators/ParticleIteratorWrapper.h"
#include "autopas/selectors/AutoTuner.h"
#include "autopas/utils/logging/Logger.h"
#include "autopas/utils/markParticleAsDeleted.h"

namespace autopas {

/**
 * The LogicHandler takes care of the containers s.t. they are all in the same valid state.
 * This is mainly done by incorporating a global container rebuild frequency, which defines when containers and their
 * neighbor lists will be rebuild.
 */
template <typename Particle>
class LogicHandler {
 public:
  /**
   * Constructor of the LogicHandler.
   * @param autoTuner
   * @param rebuildFrequency
   */
  LogicHandler(autopas::AutoTuner<Particle> &autoTuner, unsigned int rebuildFrequency)
      : _neighborListRebuildFrequency{rebuildFrequency}, _autoTuner(autoTuner) {
    checkMinimalSize();
  }

  /**
   * Collects leaving particles from buffer and potentially inserts owned particles to the container.
   * @param insertOwnedParticlesToContainer Decides whether to insert owned particles to the container.
   * @return Leaving particles.
   */
  [[nodiscard]] std::vector<Particle> collectLeavingParticlesFromBuffer(bool insertOwnedParticlesToContainer) {
    std::vector<Particle> leavingBufferParticles;
    if (insertOwnedParticlesToContainer) {
      for (auto &&p : _particleBuffer) {
        if (p.isDummy()) {
          continue;
        }
        if (utils::inBox(p.getR(), _autoTuner.getContainer()->getBoxMin(), _autoTuner.getContainer()->getBoxMax())) {
          _autoTuner.getContainer()->template addParticle(p);
        } else {
          leavingBufferParticles.push_back(p);
        }
      }
      _particleBuffer.clear();
    } else {
      for (auto iter = _particleBuffer.begin(); iter != _particleBuffer.end();) {
        auto &&p = *iter;

        auto fastRemove = [&]() {
          // Fast remove of particle, i.e., swap with last entry && pop.
          std::swap(p, leavingBufferParticles.back());
          leavingBufferParticles.pop_back();
          // Do not increment the iter afterwards!
        };
        if (p.isDummy()) {
          // We remove dummies!
          fastRemove();
        }
        if (utils::notInBox(p.getR(), _autoTuner.getContainer()->getBoxMin(), _autoTuner.getContainer()->getBoxMax())) {
          leavingBufferParticles.push_back(p);
          fastRemove();
        } else {
          ++iter;
        }
      }
    }
    return leavingBufferParticles;
  }

  /**
   * @copydoc AutoPas::updateContainer()
   */
  [[nodiscard]] std::vector<Particle> updateContainer(bool forced) {
    bool doCompleteContainerUpdate = not neighborListsAreValid() or forced;
    bool keepNeighborListsValid = not doCompleteContainerUpdate;

    if (doCompleteContainerUpdate) {
      _neighborListsAreValid = false;
    }
    // The next call also adds particles to the container if doCompleteContainerUpdate is true.
    auto leavingBufferParticles = collectLeavingParticlesFromBuffer(doCompleteContainerUpdate);

    AutoPasLog(debug, "Initiating container update.");
    auto leavingParticles = _autoTuner.getContainer()->updateContainer(keepNeighborListsValid);
    leavingParticles.insert(leavingParticles.end(), leavingBufferParticles.begin(), leavingBufferParticles.end());

    // Substract the amount of leaving particles from the number of owned particles.
    _numParticlesOwned.fetch_sub(leavingParticles.size(), std::memory_order_relaxed);
    // updateContainer deletes all halo particles.
    _haloParticleBuffer.clear();
    _numParticlesHalo.exchange(0, std::memory_order_relaxed);
    return leavingParticles;
  }

  /**
   * Pass values to the actual container.
   * @param boxMin
   * @param boxMax
   * @return Vector of particles that are outside the box after the resize.
   */
  std::vector<Particle> resizeBox(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
    const auto &oldMin = _autoTuner.getContainer()->getBoxMin();
    const auto &oldMax = _autoTuner.getContainer()->getBoxMax();

    // if nothing changed do nothing
    if (oldMin == boxMin and oldMax == boxMax) {
      return {};
    }

    // sanity check that new size is actually positive
    for (size_t i = 0; i < boxMin.size(); ++i) {
      if (boxMin[i] >= boxMax[i]) {
        utils::ExceptionHandler::exception(
            "New box size in dimension {} is not positive!\nboxMin[{}] = {}\nboxMax[{}] = {}", i, i, boxMin[i], i,
            boxMax[i]);
      }
    }

    // warn if domain changes too drastically
    const auto newLength = utils::ArrayMath::sub(boxMax, boxMin);
    const auto oldLength = utils::ArrayMath::sub(oldMax, oldMin);
    const auto relDiffLength = utils::ArrayMath::div(newLength, oldLength);
    for (size_t i = 0; i < newLength.size(); ++i) {
      // warning threshold is set arbitrary and up for change if needed
      if (relDiffLength[i] > 1.3 or relDiffLength[i] < 0.7) {
        AutoPasLog(warn,
                   "LogicHandler.resize(): Domain size changed drastically in dimension {}! Gathered AutoTuning "
                   "information might not be applicable anymore!\n"
                   "Size old box : {}\n"
                   "Size new box : {}\n"
                   "Relative diff: {}",
                   i, utils::ArrayUtils::to_string(oldLength), utils::ArrayUtils::to_string(newLength),
                   utils::ArrayUtils::to_string(relDiffLength));
      }
    }

    // check all particles
    std::vector<Particle> particlesNowOutside;
    for (auto pIter = _autoTuner.getContainer()->begin(); pIter.isValid(); ++pIter) {
      // make sure only owned ones are present
      if (not pIter->isOwned()) {
        utils::ExceptionHandler::exception(
            "LogicHandler::resizeBox() encountered non owned particle. "
            "When calling resizeBox() these should be already deleted. "
            "This could be solved by calling updateContainer() before resizeBox().");
      }
      // owned particles that are now outside are removed from the container and returned
      if (not utils::inBox(pIter->getR(), boxMin, boxMax)) {
        particlesNowOutside.push_back(*pIter);
        deleteParticle(pIter);
      }
    }

    // actually resize the container
    _autoTuner.resizeBox(boxMin, boxMax);
    _neighborListsAreValid = false;

    return particlesNowOutside;
  }

  /**
   * @copydoc AutoPas::addParticle()
   */
  void addParticle(const Particle &p) {
    if (not neighborListsAreValid()) {
      // Container has to be invalid to be able to add Particles!
      _autoTuner.getContainer()->addParticle(p);
      _numParticlesOwned.fetch_add(1, std::memory_order_relaxed);
    } else {
      // If the container is valid, we add it to the particle buffer.
      _particleBuffer.push_back(p);
    }
  }

  /**
   * @copydoc AutoPas::addOrUpdateHaloParticle()
   */
  void addHaloParticle(const Particle &haloParticle) {
    if (utils::inBox(haloParticle.getR(), _autoTuner.getContainer()->getBoxMin(),
                     _autoTuner.getContainer()->getBoxMax())) {
      utils::ExceptionHandler::exception("Trying to add a halo particle that is not OUTSIDE of the bounding box.\n" +
                                         haloParticle.toString());
    }

    auto container = _autoTuner.getContainer();
    if (not neighborListsAreValid()) {
      // If the neighbor lists are not valid, we can add the particle.
      container->template addHaloParticle</* checkInBox */ false>(haloParticle);
      _numParticlesHalo.fetch_add(1, std::memory_order_relaxed);
    } else {
      // Check if we can update an existing halo(dummy) particle.
      bool updated = _autoTuner.getContainer()->updateHaloParticle(haloParticle);
      if (not updated) {
        // If we couldn't find an existing particle, add it to the halo particle buffer.
        _haloParticleBuffer.push_back(haloParticle);
        _haloParticleBuffer.back().setOwnershipState(OwnershipState::halo);
      }
    }
  }

  /**
   * @copydoc AutoPas::deleteAllParticles()
   */
  void deleteAllParticles() {
    _neighborListsAreValid = false;
    _autoTuner.getContainer()->deleteAllParticles();
    _particleBuffer.clear();
    _haloParticleBuffer.clear();
    // all particles are gone -> reset counters.
    _numParticlesOwned.exchange(0, std::memory_order_relaxed);
    _numParticlesHalo.exchange(0, std::memory_order_relaxed);
  }

  /**
   * Deletes a single particle and updates internal particle counters.
   * @param iter
   */
  void deleteParticle(ParticleIteratorWrapper<Particle, true> &iter) {
    if ((*iter).isOwned()) {
      _numParticlesOwned.fetch_sub(1, std::memory_order_relaxed);
    } else {
      _numParticlesHalo.fetch_sub(1, std::memory_order_relaxed);
    }
    internal::markParticleAsDeleted(*iter);
  }

  /**
   * @copydoc AutoPas::iteratePairwise()
   */
  template <class Functor>
  bool iteratePairwise(Functor *f) {
    const bool doRebuild = not neighborListsAreValid();

    bool result = _autoTuner.iteratePairwise(f, doRebuild, _particleBuffer, _haloParticleBuffer);

    if (doRebuild /*we have done a rebuild now*/) {
      // list is now valid
      _neighborListsAreValid = true;
      _stepsSinceLastListRebuild = 0;
    }
    ++_stepsSinceLastListRebuild;

    return result;
  }

  /**
   * @copydoc AutoPas::begin()
   */
  autopas::ParticleIteratorWrapper<Particle, true> begin(IteratorBehavior behavior) {
    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return _autoTuner.getContainer()->begin(behavior);
  }

  /**
   * @copydoc AutoPas::begin()
   */
  autopas::ParticleIteratorWrapper<Particle, false> begin(IteratorBehavior behavior) const {
    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return std::as_const(_autoTuner).getContainer()->begin(behavior);
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  autopas::ParticleIteratorWrapper<Particle, true> getRegionIterator(std::array<double, 3> lowerCorner,
                                                                     std::array<double, 3> higherCorner,
                                                                     IteratorBehavior behavior) {
    // sanity check: Most of our stuff depends on `inBox` which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner[d] > higherCorner[d]) {
        autopas::utils::ExceptionHandler::exception(
            "Requesting region Iterator where the upper corner is lower than the lower corner!\n"
            "Lower corner: {}\n"
            "Upper corner: {}",
            lowerCorner, higherCorner);
      }
    }

    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return _autoTuner.getContainer()->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  autopas::ParticleIteratorWrapper<Particle, false> getRegionIterator(std::array<double, 3> lowerCorner,
                                                                      std::array<double, 3> higherCorner,
                                                                      IteratorBehavior behavior) const {
    // sanity check: Most of our stuff depends on `inBox` which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner[d] > higherCorner[d]) {
        autopas::utils::ExceptionHandler::exception(
            "Requesting region Iterator where the upper corner is lower than the lower corner!\n"
            "Lower corner: {}\n"
            "Upper corner: {}",
            lowerCorner, higherCorner);
      }
    }

    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return std::as_const(_autoTuner).getContainer()->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * Get the number of owned particles.
   * @return
   */
  [[nodiscard]] unsigned long getNumParticlesOwned() const { return _numParticlesOwned; }

  /**
   * Get the number of halo particles.
   * @return
   */
  [[nodiscard]] unsigned long getNumParticlesHalo() const { return _numParticlesHalo; }

 private:
  void checkMinimalSize() {
    auto container = _autoTuner.getContainer();
    // check boxSize at least cutoff + skin
    for (unsigned int dim = 0; dim < 3; ++dim) {
      if (container->getBoxMax()[dim] - container->getBoxMin()[dim] < container->getCutoff() + container->getSkin()) {
        autopas::utils::ExceptionHandler::exception(
            "Box (boxMin[{}]={} and boxMax[{}]={}) is too small.\nHas to be at least cutoff({}) + skin({}) = {}.", dim,
            container->getBoxMin()[dim], dim, container->getBoxMax()[dim], container->getCutoff(), container->getSkin(),
            container->getCutoff() + container->getSkin());
      }
    }
  }

  bool neighborListsAreValid() {
    if (_stepsSinceLastListRebuild >= _neighborListRebuildFrequency or _autoTuner.willRebuild()) {
      _neighborListsAreValid = false;
    }
    return _neighborListsAreValid;
  }

  /**
   * Specifies after how many pair-wise traversals the neighbor lists (if they exist) are to be rebuild.
   */
  const unsigned int _neighborListRebuildFrequency;

  /**
   * Reference to the AutoTuner that owns the container, ...
   */
  autopas::AutoTuner<Particle> &_autoTuner;

  /**
   * Specifies if the neighbor list is valid.
   */
  bool _neighborListsAreValid{false};

  /**
   * Steps since last rebuild
   */
  unsigned int _stepsSinceLastListRebuild{std::numeric_limits<unsigned int>::max()};

  /**
   * Atomic tracker of the number of owned particles.
   */
  std::atomic<size_t> _numParticlesOwned{0ul};

  /**
   * Atomic tracker of the number of halo particles.
   */
  std::atomic<size_t> _numParticlesHalo{0ul};

  /**
   * Buffer to store particles that should not yet be added to the container.
   * @note This buffer could potentially be replaced by a ParticleCell.
   */
  std::vector<Particle> _particleBuffer;

  /**
   * Buffer to store halo particles that should not yet be added to the container.
   * @note This buffer could potentially be replaced by a ParticleCell.
   */
  std::vector<Particle> _haloParticleBuffer;
};
}  // namespace autopas
