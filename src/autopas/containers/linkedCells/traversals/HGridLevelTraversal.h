/**
 * @file HGridLevelTraversal.h
 * @author atacann
 * @date 09.12.2024
 */

#pragma once

#include "LCTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

template <class ParticleCell>
class HGridLevelTraversal : public TraversalInterface {
 public:
  explicit HGridLevelTraversal(DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3) {}

  virtual std::unique_ptr<CellTraversal<ParticleCell>> generateNewTraversal(
      double cutoff, double interactionLength, const std::array<unsigned long, 3> &cellsPerDimension,
      const std::array<double, 3> &cellLength) = 0;

};
}  // namespace autopas
