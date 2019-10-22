/**
 * @file FmmParticle.h
 * @date 15.09.19
 * @author Joachim Marin
 */

#pragma once

#include "autopas/AutoPas.h"

class FmmParticle : public autopas::Particle {
 public:
  double charge;
  double resultFMM = 0;
  double resultExact = 0;
  double longRange = 0;
  double shortRange = 0;

  FmmParticle(std::array<double, 3> r, std::array<double, 3> v, int id, double charge) : ParticleBase(r, v, id) {
    this->charge = charge;
  }

  FmmParticle() : ParticleBase() { this->charge = 1; }
};
