//
// Created by markus on 9/25/25.
//

#pragma once
#include "PairwiseFunctorWrapperInterface.h"
#include "FunctorWrapper.h"
#include "PairwiseFunctor.h"

template <class Particle_T, class Functor_T>
class PairwiseFunctorWrapper : public FunctorWrapper<Particle_T, Functor_T>, public PairwiseFunctorWrapperInterface<Particle_T> {
public:

  PairwiseFunctorWrapper(Functor_T f) :  FunctorWrapper<Particle_T, Functor_T>(std::move(f)) {}

  virtual ~PairwiseFunctorWrapper() = default;

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {
    this->functor.AoSFunctor(i, j, newton3);
  }

  void SoAFunctorSingle(void *soaRaw, bool newton3) {
    auto &soa = *reinterpret_cast<autopas::SoAView<typename Functor_T::SoAArraysType>*>(soaRaw);
      this->functor.SoAFunctorSingle(soa, newton3);
    }

  void SoAFunctorPair (void *soaRaw1, void *soaRaw2, bool newton3) {
    auto &soa1 = *reinterpret_cast<autopas::SoAView<typename Functor_T::SoAArraysType>*>(soaRaw1);
    auto &soa2 = *reinterpret_cast<autopas::SoAView<typename Functor_T::SoAArraysType>*>(soaRaw2);
    this->functor.SoAFunctorSingle(soa1, soa2, newton3);
  }

  void SoAFunctorVerlet(void *soaRaw, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborlist, bool newton3) {
    auto &soa = *reinterpret_cast<autopas::SoAView<typename Functor_T::SoAArraysType>*>(soaRaw);
    this->functor.SoAFunctorVerlet(soa, indexFirst, neighborlist, newton3);
  }
};