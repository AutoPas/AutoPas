/**
 * @file TouchableParticle.h
 * @author seckler
 * @date 13.01.20
 */

#pragma once
#include "autopas/particles/Particle.h"

/**
 * Class that extends the default autopas::Particle with a behavior to check how often the particle was touched.
 * This class is useful for all kinds of ParticleIterator tests.
 */
class TouchableParticle : public autopas::Particle {
 public:
  /**
   * Constructor with position and id.
   * @param pos position
   * @param id id of the particle
   */
  TouchableParticle(std::array<double, 3> pos, unsigned long id)
      : autopas::Particle(pos, {0, 0, 0}, id), _numTouched(0){};

  /**
   * Constructor with position, velocity, and id.
   * @param pos position
   * @param velocity velocity
   * @param id id of the particle
   */
  TouchableParticle(std::array<double, 3> pos, std::array<double, 3> velocity, unsigned long id)
      : autopas::Particle(pos, velocity, id), _numTouched(0){};

  /**
   * Default constructor
   */
  TouchableParticle() : TouchableParticle({0., 0., 0.}, 0ul) {}

  /**
   * Touch the particle.
   * The number of times a particle was touched is saved.
   */
  void touch() { _numTouched++; }

  /**
   * Get the number that indicates how often the particle was touch.
   * @return returns how often the particle was touched.
   */
  unsigned int getNumTouched() const { return _numTouched; }

 private:
  unsigned int _numTouched;
};