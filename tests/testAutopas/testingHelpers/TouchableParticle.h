/**
 * @file TouchableParticle.h
 * @author seckler
 * @date 13.01.20
 */

#pragma once
#include <autopas/particles/Particle.h>

class TouchableParticle : public autopas::Particle {
 public:
  TouchableParticle(std::array<double, 3> pos, unsigned long id)
      : autopas::Particle(pos, {0, 0, 0}, id), _numTouched(0){};

  TouchableParticle() : TouchableParticle({0., 0., 0.}, 0ul) {}

  void touch() { _numTouched++; }

  unsigned int getNumTouched() { return _numTouched; }

 private:
  unsigned int _numTouched;
};