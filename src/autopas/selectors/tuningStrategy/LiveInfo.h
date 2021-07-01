/**
* @file LiveInfo.h
 * @author humig
 * @date 28.06.2021
*/

#pragma once

#include "autopas/utils/ArrayMath.h"
#include "autopas/containers/ParticleContainerInterface.h"

namespace autopas {

class LiveInfo {
 public:
  LiveInfo() = default;

  template<class Particle>
  void gather(const autopas::ParticleContainerInterface<Particle>& container) {
    numParticles = container.getNumParticles();
    cutoff = container.getCutoff();
    domainSize = utils::ArrayMath::sub(container.getBoxMax(), container.getBoxMin());
  }

  size_t numParticles;
  double cutoff;
  std::array<double, 3> domainSize;
};

} // namespace autopas