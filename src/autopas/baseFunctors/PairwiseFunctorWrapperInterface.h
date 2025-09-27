//
// Created by markus on 9/25/25.
//

#pragma once
#include "FunctorWrapperInterface.h"

template <typename Particle_T>
class PairwiseFunctorWrapperInterface {
public:
  virtual ~PairwiseFunctorWrapperInterface() = default;

  virtual void AoSFunctor(Particle_T &i, Particle_T &j, bool newton) = 0;

  virtual void SoAFunctorSingle(void *soa, bool newton3) = 0;
  virtual void SoAFunctorVerlet(void *soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3) = 0;
  virtual void SoAFunctorPair(void *soa1, void *soa2, bool newton3) = 0;
};