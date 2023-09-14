/**
 * @file MockTriwiseFunctor.h
 * @author muehlhaeusser
 * @date 14.09.23
 */

#pragma once

#include <gmock/gmock.h>

#include "autopas/pairwiseFunctors/TriwiseFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
//#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

template <class Particle>
class MockTriwiseFunctor : public autopas::TriwiseFunctor<Particle, MockTriwiseFunctor<Particle>> {
 public:
  MockTriwiseFunctor() : autopas::TriwiseFunctor<Particle, MockTriwiseFunctor<Particle>>(0.){};
  // virtual void AoSFunctor(Particle &i, Particle &j, bool newton3)
  MOCK_METHOD(void, AoSFunctor, (Particle &i, Particle &j, Particle &k, bool newton3), (override));

  // virtual void SoAFunctorSingle(SoAView &soa, bool newton3)
  MOCK_METHOD(void, SoAFunctorSingle, (autopas::SoAView<typename Particle::SoAArraysType> soa, bool newton3),
              (override));

  // virtual void SoAFunctorPair(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD(void, SoAFunctorPair,
              (autopas::SoAView<typename Particle::SoAArraysType> soa1,
               autopas::SoAView<typename Particle::SoAArraysType> soa2, bool newton3),
              (override));

  // virtual void SoAFunctorPair(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD(void, SoAFunctorTriple,
              (autopas::SoAView<typename Particle::SoAArraysType> soa1,
               autopas::SoAView<typename Particle::SoAArraysType> soa2,
               autopas::SoAView<typename Particle::SoAArraysType> soa3, bool newton3),
              (override));

  // virtual void SoAFunctorVerlet(SoAView &soa, const std::vector, (override)<std::vector<size_t,
  // AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
  MOCK_METHOD(void, SoAFunctorVerlet,
              (autopas::SoAView<typename Particle::SoAArraysType> soa, size_t indexFirst,
               (const std::vector<size_t, autopas::AlignedAllocator<size_t>> &), bool newton3),
              (override));

  MOCK_METHOD(void, SoALoader,
              (autopas::FullParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset));

  MOCK_METHOD(void, SoALoader,
              (autopas::ReferenceParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset));

  MOCK_METHOD(void, SoAExtractor,
              (autopas::FullParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset));

  MOCK_METHOD(void, SoAExtractor,
              (autopas::ReferenceParticleCell<Particle> & cell, autopas::SoA<typename Particle::SoAArraysType> &soa,
               size_t offset));

  // virtual bool allowsNewton3() { return true; }
  MOCK_METHOD(bool, allowsNewton3, (), (override));

  // virtual bool allowsNonNewton3() { return false; }
  MOCK_METHOD(bool, allowsNonNewton3, (), (override));

  //  bool isRelevantForTuning() { return true; }
  MOCK_METHOD(bool, isRelevantForTuning, (), (override));
};
