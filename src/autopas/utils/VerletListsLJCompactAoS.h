#pragma once

#include <vector>

#include "AlignedAllocator.h"
#include "SoA.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

template <typename ParticleType>
class VerletListsLJCompactAoS {

public:

  struct alignas (32) compactParticle {
    double posX;
    double posY;
    double posZ;
    uint32_t ownershipState;
    uint32_t typeId;
  };

  SoA<typename ParticleType::SoAArraysType> _soa;
  std::vector<compactParticle, AlignedAllocator<compactParticle>> _compactParticles;


  void resize(std::size_t n) {
    _soa.resizeArrays(n);
    _compactParticles.resize(n);
  }
};
}