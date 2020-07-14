/**
 * @file NonConstructibleParticle.h
 * @author F. Gratl
 * @date 14/07/2020
 */

#pragma once

#include "autopas/particles/Particle.h"

/**
 * A particle class without an actual constructor (only copy, etc.).
 */
class NonConstructibleParticle : public autopas::Particle {
 public:
  /**
   * Default constructor.
   */
  NonConstructibleParticle() = default;
};