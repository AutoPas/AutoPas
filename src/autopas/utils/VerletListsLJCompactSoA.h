#pragma once

#include <vector>

#include "AlignedAllocator.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

template <typename floatType>
class VerletListsLJCompactSoA {

public:

  std::vector<uint8_t, AlignedAllocator<uint8_t>> _typeIds;
  std::vector<uint8_t, AlignedAllocator<uint8_t>> _ownershipStates;
  std::vector<floatType, AlignedAllocator<floatType>> _posX;
  std::vector<floatType, AlignedAllocator<floatType>> _posY;
  std::vector<floatType, AlignedAllocator<floatType>> _posZ;
  std::vector<floatType, AlignedAllocator<floatType>> _forceX;
  std::vector<floatType, AlignedAllocator<floatType>> _forceY;
  std::vector<floatType, AlignedAllocator<floatType>> _forceZ;

  void resize(std::size_t n) {
    _typeIds.resize(n);
    _ownershipStates.resize(n);
    _posX.resize(n);
    _posY.resize(n);
    _posZ.resize(n);
    _forceX.resize(n);
    _forceY.resize(n);
    _forceZ.resize(n);
  }
};
}  // namespace autopas