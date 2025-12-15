/**
 * @file VerletClusterListsTest.h
 * @author nguyen
 * @date 21.10.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC06Traversal.h"
#include "autopas/containers/verletClusterLists/traversals/VCLC08Traversal.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "mocks/MockPairwiseFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class VerletClusterListsTest : public AutoPasTestBase {};

class CollectParticlePairsFunctor : public autopas::PairwiseFunctor<ParticleFP64, CollectParticlePairsFunctor> {
 public:
  std::vector<std::pair<ParticleFP64 *, ParticleFP64 *>> _pairs{};
  std::array<double, 3> _min;
  std::array<double, 3> _max;

  CollectParticlePairsFunctor(double cutoff, std::array<double, 3> min, std::array<double, 3> max)
      : PairwiseFunctor(cutoff), _min(min), _max(max) {}

  void initTraversal() override { _pairs.clear(); }

  void AoSFunctor(ParticleFP64 &i, ParticleFP64 &j, bool newton3) override {
    using namespace autopas::utils::ArrayMath::literals;

    auto dist = i.getR() - j.getR();
    if (autopas::utils::ArrayMath::dot(dist, dist) > getCutoff() * getCutoff() or
        not autopas::utils::inBox(i.getR(), _min, _max) or not autopas::utils::inBox(j.getR(), _min, _max))
      return;

    AUTOPAS_OPENMP(critical) {
      _pairs.emplace_back(&i, &j);
      if (newton3) _pairs.emplace_back(&j, &i);
    };
  }

  std::string getName() override { return "CollectParticlePairsFunctor"; }

  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override { return true; }
  bool allowsNonNewton3() override { return true; }

  auto getParticlePairs() { return _pairs; }
};

#if defined(AUTOPAS_USE_OPENMP)
class CollectParticlesPerThreadFunctor
    : public autopas::PairwiseFunctor<ParticleFP64, CollectParticlesPerThreadFunctor> {
 public:
  int _currentColor{};

  std::array<std::vector<std::set<ParticleFP64 *>>, 8> _particlesPerThreadPerColor;

 public:
  CollectParticlesPerThreadFunctor() : PairwiseFunctor(0) {}

  void initTraversal() override {
    for (int i = 0; i < 8; i++) {
      _particlesPerThreadPerColor[i].resize(autopas::autopas_get_max_threads());
    }
  }

  void AoSFunctor(ParticleFP64 &i, ParticleFP64 &j, bool newton3) override {
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto threadNum = autopas::autopas_get_thread_num();
    _particlesPerThreadPerColor[_currentColor][threadNum].insert(&i);
    _particlesPerThreadPerColor[_currentColor][threadNum].insert(&j);
  }

  std::string getName() override { return "CollectParticlesPerThreadFunctor"; }

  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override { return true; }
  bool allowsNonNewton3() override { return true; }

  void nextColor(int newColor) { _currentColor = newColor; }
};

class ColoringTraversalWithColorChangeNotify
    : public autopas::VCLC06Traversal<FPCell, CollectParticlesPerThreadFunctor> {
 public:
  ColoringTraversalWithColorChangeNotify(CollectParticlesPerThreadFunctor *functor, size_t clusterSize,
                                         std::function<void(int)> whenColorChanges)
      : autopas::VCLC06Traversal<FPCell, CollectParticlesPerThreadFunctor>(functor, clusterSize,
                                                                           autopas::DataLayoutOption::aos, true) {
    _whenColorChanges = std::move(whenColorChanges);
  }

  void notifyColorChange(unsigned long newColor) override { _whenColorChanges(newColor); }

 private:
  std::function<void(int)> _whenColorChanges;
};
#endif