#ifndef AUTOPAS_AUTOPAS_H
#define AUTOPAS_AUTOPAS_H

#include "autopasIncludes.h"

template <class Particle, class ParticleCell>
class AutoPas {
 public:
  /**
   * Initialize data structure, container etc...
   */
  void init() {
    //TODO:
  }
 private:
  ParticleContainer<Particle, ParticleCell> *container;
  //TODO:
};

#endif  // AUTOPAS_AUTOPAS_H
