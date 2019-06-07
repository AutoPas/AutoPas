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
  LogicHandler(autopas::AutoTuner<Particle, ParticleCell>* autoTuner, const std::array<double, 3>& boxMin,
               const std::array<double, 3>& boxMax, double cutoff, double skin, unsigned int rebuildFrequency,
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
    return std::move(_autoTuner.getContainer()->updateContainer());
  }

  /**
   * @copydoc AutoPas::addParticle()
   */
  void addParticle(Particle& p) { _autoTuner.getContainer()->addParticle(); }

  /**
   * @copydoc AutoPas::addHaloParticle()
   */
  void addHaloParticle(Particle& haloParticle) { _autoTuner.getContainer()->addHaloParticle(); }

  /**
   * @copydoc AutoPas::deleteHaloParticles()
   */
  void deleteHaloParticles() { _autoTuner.getContainer()->deleteHaloParticles(); }

  /**
   * @copydoc AutoPas::deleteAllParticles()
   */
  void deleteAllParticles() { _autoTuner.getContainer()->deleteAllParticles(); }

  /**
   * @copydoc AutoPas::iteratePairwise()
   */
  template <class Functor>
  void iteratePairwise(Functor* f) {
    _autoTuner.iteratePairwise(f);
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
    for (unsigned int dim = 0; dim < 3; ++dim) {
      if (_boxMax[dim] - _boxMin[dim] < _cutoff + _skin) {
        AutoPasLog(error, "Box (boxMin[{}]={} and boxMax[{}]={}) is too small.", dim, _boxMin[dim], dim, _boxMax[dim]);
        AutoPasLog(error, "Has to be at least cutoff({}) + skin({}) = {}.", _cutoff, _skin, _cutoff + _skin);
        autopas::utils::ExceptionHandler::exception("Box too small.");
      }
    }
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
  autopas::AutoTuner<Particle, ParticleCell>& _autoTuner;
};
}  // namespace autopas
