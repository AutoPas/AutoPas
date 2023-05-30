/**
 * @file MockFunctor.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gmock/gmock.h>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

template <class Particle>
class MockFunctor : public autopas::Functor<Particle, MockFunctor<Particle>> {
 public:
  MockFunctor() : autopas::Functor<Particle, MockFunctor<Particle>>(0.){};
  // virtual void AoSFunctor(Particle &i, Particle &j, bool newton3)
  MOCK_METHOD(void, AoSFunctor, (Particle & i, Particle &j, bool newton3), (override));

  // virtual void SoAFunctorSingle(SoAView &soa, bool newton3)
  MOCK_METHOD(void, SoAFunctorSingle, (autopas::SoAView<typename Particle::SoAArraysType> soa, bool newton3),
              (override));

  // virtual void SoAFunctorPair(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD(void, SoAFunctorPair,
              (autopas::SoAView<typename Particle::SoAArraysType> soa,
               autopas::SoAView<typename Particle::SoAArraysType> soa2, bool newton3),
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

  // virtual bool allowsMixedNewton3CallsForGlobals() { return false; }
  MOCK_METHOD(bool, allowsMixedNewton3CallsForGlobals, (), (override));

  //  bool isRelevantForTuning() { return true; }
  MOCK_METHOD(bool, isRelevantForTuning, (), (override));
};
