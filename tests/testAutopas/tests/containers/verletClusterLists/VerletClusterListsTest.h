/**
 * @file VerletClusterListsTest.h
 * @author nguyen
 * @date 21.10.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersColoringTraversal.h"
#include "autopas/particles/Particle.h"
#include "autopas/utils/WrapOpenMP.h"
#include "mocks/MockFunctor.h"
#include "mocks/MockVerletLists.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class VerletClusterListsTest : public AutoPasTestBase {};

#if defined(AUTOPAS_OPENMP)
class CollectParticlesPerThreadFunctor
    : public autopas::Functor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> {
 public:
  static int _currentColor;
#pragma omp threadprivate(_currentColor)

 public:
  std::array<std::vector<std::set<Particle *>>, 8> _particlesPerThreadPerColor;
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

  static void nextColor(int newColor) { _currentColor = newColor; }
};

int CollectParticlesPerThreadFunctor::_currentColor = 0;

class ColoringTraversalWithColorChangeNotify
    : public autopas::VerletClustersColoringTraversal<FPCell, CollectParticlesPerThreadFunctor,
                                                      autopas::DataLayoutOption::aos, true> {
 public:
  ColoringTraversalWithColorChangeNotify(CollectParticlesPerThreadFunctor *functor,
                                         std::function<void(int)> whenColorChanges)
      : autopas::VerletClustersColoringTraversal<FPCell, CollectParticlesPerThreadFunctor,
                                                 autopas::DataLayoutOption::aos, true>(functor) {
    _whenColorChanges = std::move(whenColorChanges);
  }

  void notifyColorChange(unsigned long newColor) override { _whenColorChanges(newColor); }

 private:
  std::function<void(int)> _whenColorChanges;
};
#endif