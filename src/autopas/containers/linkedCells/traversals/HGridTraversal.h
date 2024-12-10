/**
 * @file HGridTraversal.h
 * @author atacann
 * @date 09.12.2024
 */

#pragma once

#include "HGridTraversalInterface.h"
#include "HGridLevelTraversal.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

template <class ParticleCell, class Functor>
class HGridTraversal : public HGridLevelTraversal<ParticleCell>,
                       public HGridTraversalInterface {
 public:
  explicit HGridTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGridLevelTraversal<ParticleCell>(dataLayout, useNewton3) {}

  std::unique_ptr<CellTraversal<ParticleCell>> generateNewTraversal(
      double cutoff, double interactionLength, const std::array<unsigned long, 3> &cellsPerDimension,
      const std::array<double, 3> &cellLength) override {
    auto _functor = std::make_unique<Functor>(cutoff);
    auto traversal = std::make_unique<LCC08Traversal<ParticleCell, Functor>>(
        cellsPerDimension, _functor.get(), interactionLength, cellLength, this->_dataLayout, this->_useNewton3);
    return traversal;
  };

  void traverseParticles() override {}
  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_test; };
  [[nodiscard]] bool isApplicable() const override { return true;};
  void initTraversal() override {};
  void endTraversal() override {};
};
}  // namespace autopas
