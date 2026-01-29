#pragma once

#include <vector>

#include "AlignedAllocator.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

template <typename floatType>
class VerletListLJData {

  public:

  /*struct alignas(32) ReadData {
    floatType x, y, z;
    uint32_t typeId;
    uint8_t ownershipState;
    uint8_t pad[3];
  };

  struct alignas(32) WriteData {
    floatType forceX, forceY, forceZ;
    uint64_t padding;
  };*/

  /*struct alignas(32) Positions {
    std::array<floatType, 3> positionXYZ;
    uint64_t pad;
  };*/

  std::array<floatType, 4> positionsXYZ;

  /*struct alignas(32) WriteData {
    std::array<floatType, 3> forceXYZ;
    uint64_t pad;
  };*/

  std::array<floatType, 4> forceXYZ;


  std::vector<uint16_t, AlignedAllocator<uint16_t>> _typeIds;
  std::vector<uint8_t, AlignedAllocator<uint8_t>> _ownershipStates;
  std::vector<std::array<floatType, 4>, AlignedAllocator<std::array<floatType, 4>>> _positionsXYZ;
  std::vector<std::array<floatType, 4>, autopas::AlignedAllocator<std::array<floatType, 4>>> _forcesXXZ;

  void resize(std::size_t n) {
    _typeIds.resize(n);
    _ownershipStates.resize(n);
    _positionsXYZ.resize(n);
    _forcesXXZ.resize(n);
  }

};
}  // namespace autopas
