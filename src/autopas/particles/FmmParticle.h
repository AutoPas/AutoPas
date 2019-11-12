
/**
 * @file FmmParticle.h
 * @date 13.11.19
 * @author Joachim Marin
 */

#pragma once

#include <array>
#include "autopas/particles/Particle.h"

namespace autopas::fmm {

class FmmParticle : public autopas::Particle {
 public:
  double charge;
  double resultFMM = 0;
  double resultExact = 0;
  double longRange = 0;
  double shortRange = 0;

  FmmParticle(std::array<double, 3> r, std::array<double, 3> v, int id, double charge)
      : ParticleBase(r, v, id), charge(charge) {}

  FmmParticle() : ParticleBase(), charge(1) {}
};
}  // namespace autopas::fmm
