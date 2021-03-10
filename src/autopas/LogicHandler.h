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
      : _containerRebuildFrequency{rebuildFrequency}, _autoTuner(autoTuner) {
    checkMinimalSize();
  }

  /**
   * @copydoc AutoPas::updateContainer()
   */
  [[nodiscard]] std::pair<std::vector<Particle>, bool> updateContainer(bool forced) {
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

    return particlesNowOutside;
  }

  /**
   * @copydoc AutoPas::addParticle()
   */
  void addParticle(const Particle &p) {
    if (not isContainerValid()) {
      // Container has to be invalid to be able to add Particles!
      _autoTuner.getContainer()->addParticle(p);
      _numParticlesOwned.fetch_add(1, std::memory_order_relaxed);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Adding of particles not allowed while neighborlists are still valid. Please invalidate the neighborlists "
          "by calling AutoPas::updateContainer(true). Do this on EVERY AutoPas instance, i.e., on all mpi "
          "processes!");
    }
  }

  /**
   * @copydoc AutoPas::addOrUpdateHaloParticle()
   */
  void addOrUpdateHaloParticle(const Particle &haloParticle) {
    using ::autopas::utils::ArrayMath::addScalar;
    using ::autopas::utils::ArrayMath::subScalar;

    auto container = _autoTuner.getContainer();
    if (not isContainerValid()) {
      if (not utils::inBox(haloParticle.getR(), _autoTuner.getContainer()->getBoxMin(),
                           _autoTuner.getContainer()->getBoxMax())) {
        container->template addHaloParticle</* checkInBox */ false>(haloParticle);
        _numParticlesHalo.fetch_add(1, std::memory_order_relaxed);
      } else {
        utils::ExceptionHandler::exception("Trying to add a halo particle that is not OUTSIDE of the bounding box.\n" +
                                           haloParticle.toString());
      }
    } else {
      // check if the halo particle is actually a halo particle, i.e., not too far (more than skin/2) inside of the
      // domain.
      if (not utils::inBox(haloParticle.getR(), addScalar(container->getBoxMin(), container->getSkin() / 2),
                           subScalar(container->getBoxMax(), container->getSkin() / 2))) {
        bool updated = _autoTuner.getContainer()->updateHaloParticle(haloParticle);
        if (not updated) {
          // a particle has to be updated if it is within cutoff + skin/2 of the bounding box
          double dangerousDistance = container->getCutoff() + container->getSkin() / 2;

          auto dangerousBoxMin = subScalar(container->getBoxMin(), dangerousDistance);
          auto dangerousBoxMax = addScalar(container->getBoxMax(), dangerousDistance);
          bool dangerous = utils::inBox(haloParticle.getR(), dangerousBoxMin, dangerousBoxMax);
          if (dangerous) {
            // throw exception, rebuild frequency not high enough / skin too small!
            utils::ExceptionHandler::exception(
                "LogicHandler::addHaloParticle: wasn't able to update halo particle that is too close to "
                "domain (more than cutoff + skin/2). Rebuild frequency not high enough / skin too small!\n"
                "Cutoff       : {}\n" +
                "Skin         : {}\n" +
                "BoxMin       : {}\n" +
                "BoxMax       : {}\n" +
                "Dangerous Min: {}\n" +
                "Dangerous Max: {}\n" +
                "Particle     : {}\n" +
                container->getCutoff(), container->getSkin(), container->getBoxMin(), container->getBoxMax(), dangerousBoxMin, dangerousBoxMax, haloParticle.toString());
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
    internal::markParticleAsDeleted(*iter);
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
    // sanity check: Most of our stuff depends on `inBox` which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner > higherCorner) {
        autopas::utils::ExceptionHandler::exception(
            "Requesting region Iterator where upper corner is lower than lower corner!\n"
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
  autopas::ParticleIteratorWrapper<Particle, false> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const {
    // sanity check: Most of our stuff depends on `inBox` which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner > higherCorner) {
        autopas::utils::ExceptionHandler::exception(
            "Requesting region Iterator where higherCorner corner is lower than lower corner!\n"
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
  autopas::AutoTuner<Particle> &_autoTuner;

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
