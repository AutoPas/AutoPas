#pragma once

namespace autopas::containers::kokkos::traversal {

template <class Particle_T, class Functor_T>
class KokkosVCLTraversal : public autopas::TraversalInterface {
 public:
  KokkosVCLTraversal( Functor_T* f,DataLayoutOption dataLayout, bool newton3)
      : TraversalInterface(dataLayout, newton3), _functor(f) {}

  ~KokkosVCLTraversal() override = default;

  [[nodiscard]] autopas::TraversalOption getTraversalType() const override { return autopas::TraversalOption::kk_vcl; }

  [[nodiscard]] bool isApplicable() const override {}
  void initTraversal() override {}
  void traverseParticles() override {}
  void endTraversal() override {}

 private:
  Functor_T _functor;
};

}  // namespace autopas::containers::kokkos::traversal
