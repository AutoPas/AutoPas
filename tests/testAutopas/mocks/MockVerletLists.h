/**
 * @file MockVerletLists.h
 * @author seckler
 * @date 26.04.18
 */

#pragma once

#include <gmock/gmock.h>
#include "containers/VerletLists.h"

// gmock does not write overrides, so we suppress that warning here!
#if __GNUC__ >= 5
// Disable GCC 5's -Wsuggest-override warnings in gtest
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"
#endif

template <typename Particle, typename ParticleCell>
class MockVerletLists : public autopas::VerletLists<Particle, ParticleCell> {
 public:
  MockVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin,
                  unsigned int rebuildFrequency = 1)
      : autopas::VerletLists<Particle, ParticleCell>(boxMin, boxMax, cutoff, skin, rebuildFrequency) {}
  // MOCK_METHOD2_T(iteratePairwiseAoS, void(Functor, bool));  // we are not
  // allowed to mock this! We also don't want to! as this is exactly the
  // function we want to test
  MOCK_METHOD1_T(addParticle, void(Particle& p));

  MOCK_METHOD1_T(addHaloParticle, void(Particle& haloParticle));

  MOCK_METHOD0(updateContainer, void());

 protected:
  MOCK_METHOD1(updateVerletListsAoS, void(bool));

  void addParticleVerletLists(Particle& p) { autopas::VerletLists<Particle, ParticleCell>::addParticle(p); }
  void addHaloParticleVerletLists(Particle& p) { autopas::VerletLists<Particle, ParticleCell>::addHaloParticle(p); }
  void updateContainerVerletLists() { autopas::VerletLists<Particle, ParticleCell>::updateContainer(); }

  friend class VerletListsTest_testRebuildFrequencyAlways_Test;
  friend class VerletListsTest_testRebuildFrequencyEvery3_Test;
  friend class VerletListsTest_testForceRebuild_Test;
};

#if __GNUC__ >= 5
#pragma GCC diagnostic pop
#endif