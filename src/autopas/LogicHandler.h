/**
 * @file LogicHandler.h
 * @author seckler
 * @date 31.05.19
 */

#pragma once
#include <limits>
#include "autopas/iterators/ParticleIteratorWrapper.h"
#include "autopas/selectors/AutoTuner.h"
#include "autopas/utils/Logger.h"

namespace autopas {

/**
 * The LogicHandler takes care of the containers s.t. they are all in the same valid state.
 * This is mainly done by incorporating a global container rebuild frequency, which defines when containers and their
 * neighbor lists will be rebuild.
 */
template <typename Particle, typename ParticleCell>
class LogicHandler {
 public:
  /**
   * Constructor of the LogicHandler.
   * @param autoTuner
   * @param rebuildFrequency
   */
  LogicHandler(autopas::AutoTuner<Particle, ParticleCell> &autoTuner, unsigned int rebuildFrequency)
      : _containerRebuildFrequency{rebuildFrequency}, _autoTuner(autoTuner) {
    checkMinimalSize();
  }

  /**
   * @copydoc AutoPas::updateContainer()
   * @param forced specifies whether an update of the container is enforced.
   */
  AUTOPAS_WARN_UNUSED_RESULT
  std::pair<std::vector<Particle>, bool> updateContainer(bool forced) {
    if (not isContainerValid() or forced) {
      AutoPasLog(debug, "Initiating container update.");
      _containerIsValid = false;
      auto returnPair = std::make_pair(std::move(_autoTuner.getContainer()->updateContainer()), true);
      // update container returns the particles which were previously owned and are now removed.
      // Therefore remove them from the counter.
      _numParticlesOwned.fetch_sub(returnPair.first.size(), std::memory_order_relaxed);
      // updateContainer deletes all halo particles.
      _numParticlesHalo.exchange(0, std::memory_order_relaxed);
      return returnPair;
    } else {
      AutoPasLog(debug, "Skipping container update.");
      return std::make_pair(std::vector<Particle>{}, false);
    }
  }

  /**
   * @copydoc AutoPas::addParticle()
   */
  void addParticle(Particle &p) {
    if (not isContainerValid()) {
      _autoTuner.getContainer()->addParticle(p);
      _numParticlesOwned.fetch_add(1, std::memory_order_relaxed);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Adding of particles not allowed while neighborlists are still valid. Please invalidate the neighborlists "
          "by calling AutoPas::updateContainerForced(). Do this on EVERY AutoPas instance, i.e., on all mpi "
          "processes!");
    }
  }

  /**
   * @copydoc AutoPas::addOrUpdateHaloParticle()
   */
  void addOrUpdateHaloParticle(Particle &haloParticle) {
    auto container = _autoTuner.getContainer();
    if (not isContainerValid()) {
      if (not utils::inBox(haloParticle.getR(), _autoTuner.getContainer()->getBoxMin(),
                           _autoTuner.getContainer()->getBoxMax())) {
        container->addHaloParticle(haloParticle);
        _numParticlesHalo.fetch_add(1, std::memory_order_relaxed);
      } else {
        utils::ExceptionHandler::exception("Trying to add a halo particle that is not OUTSIDE of the bounding box.\n" +
                                           haloParticle.toString());
      }
    } else {
      // check if the halo particle is actually a halo particle, i.e., not too far (more than skin/2) inside of the
      // domain.
      if (not utils::inBox(haloParticle.getR(),
                           utils::ArrayMath::addScalar(container->getBoxMin(), container->getSkin() / 2),
                           utils::ArrayMath::subScalar(container->getBoxMax(), container->getSkin() / 2))) {
        bool updated = _autoTuner.getContainer()->updateHaloParticle(haloParticle);
        if (not updated) {
          // a particle has to be updated if it is within cutoff + skin/2 of the bounding box
          double dangerousDistance = container->getCutoff() + container->getSkin() / 2;

          bool dangerous =
              utils::inBox(haloParticle.getR(), utils::ArrayMath::subScalar(container->getBoxMin(), dangerousDistance),
                           utils::ArrayMath::addScalar(container->getBoxMax(), dangerousDistance));
          if (dangerous) {
            // throw exception, rebuild frequency not high enough / skin too small!
            utils::ExceptionHandler::exception(
                "LogicHandler::addHaloParticle: wasn't able to update halo particle that is too close to "
                "domain (more than cutoff + skin/2). Rebuild frequency not high enough / skin too small!");
          }
        }
      } else {
        // throw exception, rebuild frequency not high enough / skin too small!
        utils::ExceptionHandler::exception(
            "LogicHandler::addHaloParticle: trying to update halo particle that is too far inside domain "
            "(more than skin/2). Rebuild frequency not high enough / skin too small for particle \n" +
            haloParticle.toString());
      }
    }
  }

  /**
   * @copydoc AutoPas::deleteAllParticles()
   */
  void deleteAllParticles() {
    _containerIsValid = false;
    _autoTuner.getContainer()->deleteAllParticles();
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
    iter.deleteCurrentParticle();
  }

  /**
   * @copydoc AutoPas::iteratePairwise()
   */
  template <class Functor>
  bool iteratePairwise(Functor *f) {
    const bool doRebuild = not isContainerValid();
    bool result = _autoTuner.iteratePairwise(f, doRebuild);
    if (doRebuild /*we have done a rebuild now*/) {
      // list is now valid
      _containerIsValid = true;
      _stepsSinceLastContainerRebuild = 0;
    }
    ++_stepsSinceLastContainerRebuild;

    return result;
  }

  /**
   * @copydoc AutoPas::begin()
   */
  autopas::ParticleIteratorWrapper<Particle, true> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return _autoTuner.getContainer()->begin(behavior);
  }

  /**
   * @copydoc AutoPas::begin()
   */
  autopas::ParticleIteratorWrapper<Particle, false> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const {
    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return std::as_const(_autoTuner).getContainer()->begin(behavior);
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  autopas::ParticleIteratorWrapper<Particle, true> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return _autoTuner.getContainer()->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  autopas::ParticleIteratorWrapper<Particle, false> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const {
    /// @todo: we might have to add a rebuild here, if the verlet cluster lists are used.
    return _autoTuner.getContainer()->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * Get the number of owned particles.
   * @return
   */
  unsigned long getNumParticlesOwned() const { return _numParticlesOwned; }

  /**
   * Get the number of halo particles.
   * @return
   */
  unsigned long getNumParticlesHalo() const { return _numParticlesHalo; }

 private:
  void checkMinimalSize() {
    auto container = _autoTuner.getContainer();
    // check boxSize at least cutoff + skin
    for (unsigned int dim = 0; dim < 3; ++dim) {
      if (container->getBoxMax()[dim] - container->getBoxMin()[dim] < container->getCutoff() + container->getSkin()) {
        AutoPasLog(error, "Box (boxMin[{}]={} and boxMax[{}]={}) is too small.", dim, container->getBoxMin()[dim], dim,
                   container->getBoxMax()[dim]);
        AutoPasLog(error, "Has to be at least cutoff({}) + skin({}) = {}.", container->getCutoff(),
                   container->getSkin(), container->getCutoff() + container->getSkin());
        autopas::utils::ExceptionHandler::exception("Box too small.");
      }
    }
  }

  bool isContainerValid() {
    if (_stepsSinceLastContainerRebuild >= _containerRebuildFrequency or _autoTuner.willRebuild()) {
      _containerIsValid = false;
    }
    return _containerIsValid;
  }

  /**
   * Specifies after how many pair-wise traversals the container and their neighbor lists (if they exist) are to be
   * rebuild.
   */
  const unsigned int _containerRebuildFrequency;

  /**
   * Reference to the AutoTuner that owns the container, ...
   */
  autopas::AutoTuner<Particle, ParticleCell> &_autoTuner;

  /**
   * Specifies if the neighbor list is valid.
   */
  bool _containerIsValid{false};

  /**
   * Steps since last rebuild
   */
  unsigned int _stepsSinceLastContainerRebuild{std::numeric_limits<unsigned int>::max()};

  /**
   * Atomic tracker of the number of owned particles.
   */
  std::atomic<size_t> _numParticlesOwned{0ul};

  /**
   * Atomic tracker of the number of halo particles.
   */
  std::atomic<size_t> _numParticlesHalo{0ul};
};
}  // namespace autopas
