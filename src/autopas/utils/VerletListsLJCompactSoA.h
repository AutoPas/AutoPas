#pragma once

#include <vector>

#include "AlignedAllocator.h"
#include "SoA.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

template <typename ParticleType>
class VerletListsLJCompactSoA {

public:

  SoA<typename ParticleType::SoAArraysType> _soa;
  std::vector<uint8_t, AlignedAllocator<uint8_t>> _typeIds;
  std::vector<uint8_t, AlignedAllocator<uint8_t>> _ownershipStates;

  void resize(std::size_t n) {
    _soa.resizeArrays(n);
    _typeIds.resize(n);
    _ownershipStates.resize(n);
  }
};
}