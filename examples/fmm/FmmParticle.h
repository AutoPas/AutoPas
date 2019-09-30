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
  double resultFMM;
  double resultExact;

  FmmParticle(std::array<double, 3> r, std::array<double, 3> v, int id, double charge) : ParticleBase(r, v, id) {
    this->charge = charge;
    this->resultFMM = 0;
    this->resultExact = 0;
  }

  FmmParticle() : ParticleBase() {
    this->charge = 1;
    this->resultFMM = 0;
    this->resultExact = 0;
  }
};
