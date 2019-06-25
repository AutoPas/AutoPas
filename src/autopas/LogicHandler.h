/**
 * @file LogicHandler.h
 * @author seckler
 * @date 31.05.19
 */

#pragma once
#include "autopas/selectors/AutoTuner.h"
#include "autopas/utils/Logger.h"

namespace autopas {

/**
 * The LogicHandler takes care of the containers s.t. they are all in the same valid state.
 */
template <typename Particle, typename ParticleCell>
class LogicHandler {
 public:
  /**
   * Constructor of the LogicHandler.
   * @param autoTuner
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param rebuildFrequency
   * @param tuningInterval
   * @param numSamples
   */
  LogicHandler(autopas::AutoTuner<Particle, ParticleCell> *autoTuner, const std::array<double, 3> &boxMin,
               const std::array<double, 3> &boxMax, double cutoff, double skin, unsigned int rebuildFrequency,
               unsigned int tuningInterval, unsigned int numSamples)
      : _boxMin{boxMin},
        _boxMax{boxMax},
        _cutoff{cutoff},
        _skin{skin},
        _verletRebuildFrequency{rebuildFrequency},
        _tuningInterval{tuningInterval},
        _numSamples{numSamples},
        _autoTuner(*autoTuner) {
    doAssertions();
  }

  /**
   * @copydoc AutoPas::updateContainer()
   */
  std::vector<Particle> AUTOPAS_WARN_UNUSED_RESULT updateContainer() {
    if (not isNeighborListValid()) {
      AutoPasLog(debug, "Initiating container update.");
      _neighborListIsValid = false;
      return std::move(_autoTuner.getContainer()->updateContainer());
    } else {
      AutoPasLog(debug, "Skipping container update.");
      return std::vector<Particle>{};
    }
  }

  /**
   * @copydoc AutoPas::addParticle()
   */
  void addParticle(Particle &p) {
    if (not isNeighborListValid()) {
      _autoTuner.getContainer()->addParticle(p);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Adding of particles not allowed while neighborlists are still valid. Please invalidate the neighborlists "
          "by calling AutoPas::invalidateLists(). Do this on EVERY AutoPas instance, i.e., on all mpi processes!");
    }
  }

  /**
   * @copydoc AutoPas::addHaloParticle()
   */
  void addOrUpdateHaloParticle(Particle &haloParticle) {
    if (not isNeighborListValid()) {
      _autoTuner.getContainer()->addHaloParticle(haloParticle);
    } else {
      if (not utils::inBox(haloParticle.getR(), ArrayMath::addScalar(_boxMin, _skin / 2),
                           ArrayMath::subScalar(_boxMax, _skin / 2))) {
        bool updated = _autoTuner.getContainer()->updateHaloParticle(haloParticle);
        if (not updated) {
          // a particle has to be updated if it is within cutoff + skin/2 of the bounding box
          double dangerousDistance = _cutoff + _skin / 2;

          bool dangerous = utils::inBox(haloParticle.getR(), ArrayMath::subScalar(_boxMin, dangerousDistance),
                                        ArrayMath::addScalar(_boxMax, dangerousDistance));
          if (dangerous) {
            // throw exception, rebuild frequency not high enough / skin too small!
            utils::ExceptionHandler::exception(
                "VerletListsLinkedBase::addHaloParticle: wasn't able to update halo particle that is too close to "
                "domain (more than cutoff + skin/2). Rebuild frequency not high enough / skin too small!");
          }
        }
      } else {
        // throw exception, rebuild frequency not high enough / skin too small!
        utils::ExceptionHandler::exception(
            "VerletListsLinkedBase::addHaloParticle: trying to update halo particle that is too far inside domain "
            "(more than skin/2). Rebuild frequency not high enough / skin too small!");
      }
    }
  }

  /**
   * @copydoc AutoPas::deleteHaloParticles()
   */
  void deleteHaloParticles() {
    _neighborListIsValid = false;
    _autoTuner.getContainer()->deleteHaloParticles();
  }

  /**
   * @copydoc AutoPas::deleteAllParticles()
   */
  void deleteAllParticles() {
    _neighborListIsValid = false;
    _autoTuner.getContainer()->deleteAllParticles();
  }

  /**
   * @copydoc AutoPas::iteratePairwise()
   */
  template <class Functor>
  void iteratePairwise(Functor *f) {
    const bool doRebuild = not isNeighborListValid();
    _autoTuner.iteratePairwise(f, doRebuild);
    if (doRebuild /*we have done a rebuild now*/) {
      // list is now valid
      _neighborListIsValid = true;
      _stepsSinceLastContainerRebuild = 0;
    }
    ++_stepsSinceLastContainerRebuild;
  }

  /**
   * @copydoc AutoPas::begin()
   */
  autopas::ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    return _autoTuner.getContainer()->begin(behavior);
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  autopas::ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    return _autoTuner.getContainer()->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

 private:
  void doAssertions() {
    // check boxSize at least cutoff + skin
    for (unsigned int dim = 0; dim < 3; ++dim) {
      if (_boxMax[dim] - _boxMin[dim] < _cutoff + _skin) {
        AutoPasLog(error, "Box (boxMin[{}]={} and boxMax[{}]={}) is too small.", dim, _boxMin[dim], dim, _boxMax[dim]);
        AutoPasLog(error, "Has to be at least cutoff({}) + skin({}) = {}.", _cutoff, _skin, _cutoff + _skin);
        autopas::utils::ExceptionHandler::exception("Box too small.");
      }
    }
  }

  bool isNeighborListValid() {
    return _neighborListIsValid and _stepsSinceLastContainerRebuild < _verletRebuildFrequency and
           not _autoTuner.willRebuild();
  }

  /**
   * Lower corner of the container.
   */
  std::array<double, 3> _boxMin;

  /**
   * Upper corner of the container.
   */
  std::array<double, 3> _boxMax;

  /**
   * Cutoff radius to be used in this container.
   */
  double _cutoff;

  /**
   * Skin.
   */
  double _skin;

  /**
   * Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   */
  unsigned int _verletRebuildFrequency;

  /**
   * Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  unsigned int _tuningInterval;

  /**
   * Number of samples the tuner should collect for each combination.
   */
  unsigned int _numSamples;

  /**
   * Reference to the AutoTuner that owns the container, ...
   */
  autopas::AutoTuner<Particle, ParticleCell> &_autoTuner;

  /**
   * Specifies if the neighbor list is valid.
   */
  bool _neighborListIsValid{false};

  /**
   * Steps since last rebuild
   */
  unsigned int _stepsSinceLastContainerRebuild{UINT32_MAX};
};
}  // namespace autopas
