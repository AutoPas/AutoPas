/**
 * @file VerletClusterListsTest.h
 * @author nguyen
 * @date 21.10.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersColoringTraversal.h"
#include "autopas/particles/Particle.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class VerletClusterListsTest : public AutoPasTestBase {};

class CollectParticlePairsFunctor
    : public autopas::Functor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> {
 public:
  std::vector<std::pair<Particle *, Particle *>> _pairs{};
  std::array<double, 3> _min;
  std::array<double, 3> _max;

 public:
  CollectParticlePairsFunctor(double cutoff, std::array<double, 3> min, std::array<double, 3> max)
      : Functor(cutoff), _min(min), _max(max) {}

  void initTraversal() override { _pairs.clear(); }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    auto dist = autopas::utils::ArrayMath::sub(i.getR(), j.getR());
    if (autopas::utils::ArrayMath::dot(dist, dist) > getCutoff() * getCutoff() or
        not autopas::utils::inBox(i.getR(), _min, _max) or not autopas::utils::inBox(j.getR(), _min, _max))
      return;

#if defined(AUTOPAS_OPENMP)
#pragma omp critical
#endif
    {
      _pairs.emplace_back(&i, &j);
      if (newton3) _pairs.emplace_back(&j, &i);
    };
  }

  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override { return true; }
  bool allowsNonNewton3() override { return true; }

  bool isAppropriateClusterSize(unsigned int clusterSize, autopas::DataLayoutOption::Value dataLayout) const override {
    return true;
  }

  auto getParticlePairs() { return _pairs; }
};

#if defined(AUTOPAS_OPENMP)
class CollectParticlesPerThreadFunctor
    : public autopas::Functor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> {
 public:
  static int _currentColor;
#pragma omp threadprivate(_currentColor)

  std::array<std::vector<std::set<Particle *>>, 8> _particlesPerThreadPerColor;

 public:
  CollectParticlesPerThreadFunctor() : Functor(0) {}

  void initTraversal() override {
    for (int i = 0; i < 8; i++) {
      _particlesPerThreadPerColor[i].resize(autopas::autopas_get_max_threads());
    }
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    auto threadNum = autopas::autopas_get_thread_num();
    _particlesPerThreadPerColor[_currentColor][threadNum].insert(&i);
    _particlesPerThreadPerColor[_currentColor][threadNum].insert(&j);
  }

  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override { return true; }
  bool allowsNonNewton3() override { return true; }

  bool isAppropriateClusterSize(unsigned int clusterSize, autopas::DataLayoutOption::Value dataLayout) const override {
    return dataLayout == autopas::DataLayoutOption::aos;  // this functor supports clusters only for aos!
  }

  static void nextColor(int newColor) { _currentColor = newColor; }
};

int CollectParticlesPerThreadFunctor::_currentColor = 0;

class ColoringTraversalWithColorChangeNotify
    : public autopas::VerletClustersColoringTraversal<FPCell, CollectParticlesPerThreadFunctor,
                                                      autopas::DataLayoutOption::aos, true> {

 public:
  ColoringTraversalWithColorChangeNotify(CollectParticlesPerThreadFunctor *functor, size_t clusterSize,
                                         std::function<void(int)> whenColorChanges)
      : autopas::VerletClustersColoringTraversal<FPCell, CollectParticlesPerThreadFunctor,
                                                 autopas::DataLayoutOption::aos, true>(functor, clusterSize) {
    _whenColorChanges = std::move(whenColorChanges);
  }

  void notifyColorChange(unsigned long newColor) override { _whenColorChanges(newColor); }

 private:
  std::function<void(int)> _whenColorChanges;
};
#endif