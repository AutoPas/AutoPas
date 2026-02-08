/**
 * @file MockPairwiseFunctor.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gmock/gmock.h>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

template <class Particle_T>
class MockPairwiseFunctor : public autopas::PairwiseFunctor<Particle_T, MockPairwiseFunctor<Particle_T>> {
 public:
  MockPairwiseFunctor() : autopas::PairwiseFunctor<Particle_T, MockPairwiseFunctor<Particle_T>>(0.){};
  // virtual void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3)
  MOCK_METHOD(void, AoSFunctor, (Particle_T & i, Particle_T &j, bool newton3), (override));

  // virtual void SoAFunctorSingle(SoAView &soa, bool newton3)
  MOCK_METHOD(void, SoAFunctorSingle, (autopas::SoAView<typename Particle_T::SoAArraysType> soa, bool newton3),
              (override));

  // virtual void SoAFunctorPair(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD(void, SoAFunctorPair,
              (autopas::SoAView<typename Particle_T::SoAArraysType> soa,
               autopas::SoAView<typename Particle_T::SoAArraysType> soa2, bool newton3),
              (override));

  // virtual void SoAFunctorVerlet(SoAView &soa, const std::vector, (override)<std::vector<size_t,
  // AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
  MOCK_METHOD(void, SoAFunctorVerlet,
              (autopas::SoAView<typename Particle_T::SoAArraysType> soa, size_t indexFirst,
               (const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &), bool newton3),
              (override));

  MOCK_METHOD(void, SoAFunctorVerletOptimized,
            (autopas::SoAView<typename Particle_T::SoAArraysType> soa, size_t indexFirst,
             (const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &), bool newton3),
            (override));

  MOCK_METHOD(void, SoALoader,
              (autopas::FullParticleCell<Particle_T> & cell, autopas::SoA<typename Particle_T::SoAArraysType> &soa,
               size_t offset, bool skipSoAResize));

  MOCK_METHOD(void, SoALoader,
              (autopas::ReferenceParticleCell<Particle_T> & cell, autopas::SoA<typename Particle_T::SoAArraysType> &soa,
               size_t offset, bool skipSoAResize));

  MOCK_METHOD(void, SoAExtractor,
              (autopas::FullParticleCell<Particle_T> & cell, autopas::SoA<typename Particle_T::SoAArraysType> &soa,
               size_t offset));

  MOCK_METHOD(void, SoAExtractor,
              (autopas::ReferenceParticleCell<Particle_T> & cell, autopas::SoA<typename Particle_T::SoAArraysType> &soa,
               size_t offset));

  // virtual bool allowsNewton3() { return true; }
  MOCK_METHOD(bool, allowsNewton3, (), (override));

  // virtual bool allowsNonNewton3() { return false; }
  MOCK_METHOD(bool, allowsNonNewton3, (), (override));

  //  bool isRelevantForTuning() { return true; }
  MOCK_METHOD(bool, isRelevantForTuning, (), (override));

  //  std::string getName() { return "functorName"; }
  MOCK_METHOD(std::string, getName, (), (override));
};
