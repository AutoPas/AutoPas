/**
 * @file MockFunctor.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gmock/gmock.h>
#include "autopasIncludes.h"

// gmock does not write overrides, so we suppress that warning here!
#if __GNUC__ >= 5
// Disable GCC 5's -Wsuggest-override warnings in gtest
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"
#endif

template <class Particle, class ParticleCell>
class MockFunctor : public autopas::Functor<Particle, ParticleCell> {
 public:
  // virtual void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) {}
  MOCK_METHOD2_T(AoSFunctor, void(Particle &i, Particle &j));
  MOCK_METHOD3_T(AoSFunctor, void(Particle &i, Particle &j, bool newton3));

  // virtual void SoAFunctor(SoA &soa, bool newton3 = true) {}
  MOCK_METHOD1(SoAFunctor, void(autopas::SoA &soa));
  MOCK_METHOD2(SoAFunctor, void(autopas::SoA &soa, bool newton3));

  // virtual void SoAFunctor(SoA &soa1, SoA &soa2, bool newton3 = true) {}
  MOCK_METHOD2(SoAFunctor, void(autopas::SoA &soa, autopas::SoA &soa2));
  MOCK_METHOD3(SoAFunctor,
               void(autopas::SoA &soa, autopas::SoA &soa2, bool newton3));

  // virtual void SoAFunctor(SoA &soa, std::vector<std::vector<size_t>>
  // &neighborList, size_t iFrom, size_t iTo, bool newton3 = true{})
  MOCK_METHOD4(SoAFunctor,
               void(autopas::SoA &soa, std::vector<std::vector<size_t>> &,
                    size_t, size_t));
  MOCK_METHOD5(SoAFunctor,
               void(autopas::SoA &soa, std::vector<std::vector<size_t>> &,
                    size_t, size_t, bool));

  // virtual void SoALoader(ParticleCell &cell, autopas::SoA *soa, size_t
  // offset=0) {}
  MOCK_METHOD2_T(SoALoader, void(ParticleCell &cell, autopas::SoA *soa));
  MOCK_METHOD3_T(SoALoader,
                 void(ParticleCell &cell, autopas::SoA *soa, size_t offset));

  // virtual void SoAExtractor(ParticleCell *cell, autopas::SoA *soa, size_t
  // offset=0) {}
  MOCK_METHOD2_T(SoAExtractor, void(ParticleCell *cell, autopas::SoA *soa));
  MOCK_METHOD3_T(SoAExtractor,
                 void(ParticleCell *cell, autopas::SoA *soa, size_t offset));

  // virtual bool allowsNewton3() { return true; }
  MOCK_METHOD0(allowsNewton3, bool());

  // virtual bool allowsNonNewton3() { return false; }
  MOCK_METHOD0(allowsNonNewton3, bool());
};

#if __GNUC__ >= 5
#pragma GCC diagnostic pop
#endif