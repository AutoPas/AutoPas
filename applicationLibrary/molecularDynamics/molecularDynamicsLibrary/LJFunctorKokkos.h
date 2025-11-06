#pragma once
#include "autopas/baseFunctors/KokkosFunctor.h"
namespace mdLib {

template <typename Particle_T>
class LJFunctorKokkos final : public autopas::KokkosFunctor<Particle_T, LJFunctorKokkos<Particle_T>> {
 public:
  template <typename Space>
  using SoAArraysType = autopas::KokkosFunctor<Particle_T, LJFunctorKokkos>::template SoAArraysType<Space>;

  explicit LJFunctorKokkos(double cutoff) : autopas::KokkosFunctor<Particle_T, LJFunctorKokkos>(cutoff) {}

  bool allowsNewton3() { return true; }
  bool allowsNonNewton3() { return false; }
  bool isRelevantForTuning() { return true; }
  std::string getName() { return "Lennard-Jones Kokkos"; }

  template <typename Space>
  KOKKOS_FUNCTION static void KokkosSoAFunctor(SoAArraysType<Space> soa, size_t index) {};
};

}  // namespace mdLib