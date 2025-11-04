#pragma once
#include "autopas/baseFunctors/PairwiseFunctor.h"

namespace mdLib {

template <typename Particle_T>
class LJFunctorKokkos : public autopas::PairwiseFunctor<Particle_T, LJFunctorKokkos<Particle_T>> {
 public:
  explicit LJFunctorKokkos(double cutoff) : autopas::PairwiseFunctor<Particle_T, LJFunctorKokkos<Particle_T>>(cutoff) {}

  void SoAFunctorSingle(autopas::SoAView<typename Particle_T::SoAArraysType> soa, bool newton3) override {}

  void SoAFunctorVerlet(autopas::SoAView<typename Particle_T::SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<unsigned long>> &neighborList,
                        bool newton3) override {}

  void SoAFunctorPair(autopas::SoAView<typename Particle_T::SoAArraysType> soa1,
                      autopas::SoAView<typename Particle_T::SoAArraysType> soa2, bool newton3) override {}

  bool allowsNewton3() { return true; }
  bool allowsNonNewton3() { return false; }
  bool isRelevantForTuning() { return true; }
  std::string getName() { return "Lennard-Jones Kokkos"; }
};

}  // namespace mdLib