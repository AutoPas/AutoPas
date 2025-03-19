/**
 * @file MockTriwiseFunctor.h
 * @author muehlhaeusser
 * @date 14.09.23
 */

#pragma once

#include <gmock/gmock.h>

#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/options/DataLayoutOption.h"

template <class ParticleT>
class MockTriwiseFunctor : public autopas::TriwiseFunctor<ParticleT, MockTriwiseFunctor<ParticleT>> {
 public:
  MockTriwiseFunctor() : autopas::TriwiseFunctor<ParticleT, MockTriwiseFunctor<ParticleT>>(0.){};
  // virtual void AoSFunctor(ParticleT &i, ParticleT &j, bool newton3)
  MOCK_METHOD(void, AoSFunctor, (ParticleT & i, ParticleT &j, ParticleT &k, bool newton3), (override));

  // virtual void SoAFunctorSingle(SoAView &soa, bool newton3)
  MOCK_METHOD(void, SoAFunctorSingle, (autopas::SoAView<typename ParticleT::SoAArraysType> soa, bool newton3),
              (override));

  // virtual void SoAFunctorPair(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD(void, SoAFunctorPair,
              (autopas::SoAView<typename ParticleT::SoAArraysType> soa1,
               autopas::SoAView<typename ParticleT::SoAArraysType> soa2, bool newton3),
              (override));

  // virtual void SoAFunctorPair(SoAView &soa1, SoAView &soa2, bool newton3)
  MOCK_METHOD(void, SoAFunctorTriple,
              (autopas::SoAView<typename ParticleT::SoAArraysType> soa1,
               autopas::SoAView<typename ParticleT::SoAArraysType> soa2,
               autopas::SoAView<typename ParticleT::SoAArraysType> soa3, bool newton3),
              (override));

  // virtual void SoAFunctorVerlet(SoAView &soa, const std::vector, (override)<std::vector<size_t,
  // AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
  MOCK_METHOD(void, SoAFunctorVerlet,
              (autopas::SoAView<typename ParticleT::SoAArraysType> soa, size_t indexFirst,
               (const std::vector<size_t, autopas::AlignedAllocator<size_t>> &), bool newton3),
              (override));

  MOCK_METHOD(void, SoALoader,
              (autopas::FullParticleCell<ParticleT> & cell, autopas::SoA<typename ParticleT::SoAArraysType> &soa,
               size_t offset, bool skipSoAResize));

  MOCK_METHOD(void, SoALoader,
              (autopas::ReferenceParticleCell<ParticleT> & cell, autopas::SoA<typename ParticleT::SoAArraysType> &soa,
               size_t offset, bool skipSoAResize));

  MOCK_METHOD(void, SoAExtractor,
              (autopas::FullParticleCell<ParticleT> & cell, autopas::SoA<typename ParticleT::SoAArraysType> &soa,
               size_t offset));

  MOCK_METHOD(void, SoAExtractor,
              (autopas::ReferenceParticleCell<ParticleT> & cell, autopas::SoA<typename ParticleT::SoAArraysType> &soa,
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
