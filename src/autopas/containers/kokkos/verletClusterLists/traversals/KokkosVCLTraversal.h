#pragma once

#include "autopas/containers/TraversalInterface.h"

namespace autopas::containers::kokkos::traversal {

template <class Particle_T, typename Functor_T>
class KokkosVCLTraversal : public autopas::TraversalInterface {
 public:
  KokkosVCLTraversal( Functor_T* f, DataLayoutOption dataLayout, bool newton3, double cutoff)
      : TraversalInterface(dataLayout, newton3), _functor(cutoff) {}

  ~KokkosVCLTraversal() override = default;

  [[nodiscard]] autopas::TraversalOption getTraversalType() const override { return autopas::TraversalOption::kk_vcl; }

  [[nodiscard]] bool isApplicable() const override {
    // TODO
    return true;
  }
  void initTraversal() override {}
  void traverseParticles() override {}
  void endTraversal() override {}

 private:
  Functor_T _functor;
};

}  // namespace autopas::containers::kokkos::traversal
