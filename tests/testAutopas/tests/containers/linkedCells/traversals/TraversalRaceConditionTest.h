/**
 * @file TraversalRaceConditionTest.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include <vector>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/containers/linkedCells/traversals/SlicedTraversal.h"
#include "autopas/utils/ArrayMath.h"
#include "testingHelpers/commonTypedefs.h"

class TraversalRaceConditionTest : public AutoPasTestBase {
 public:
  TraversalRaceConditionTest() = default;

  ~TraversalRaceConditionTest() override = default;

  void fillWithParticles(autopas::AutoPas<Particle, FPCell> &autoPas, std::array<size_t, 3> particlesPerDim);

  /*
   * Simple AoS only functor which repulses paritcles from each other with a
   * constant force of 1.
   */
  class SimpleFunctor : public autopas::Functor<Particle, FPCell> {
   public:
    using SoAArraysType = Particle::SoAArraysType;
    using ParticleCell = FPCell;
    using floatType = typename Particle::ParticleFloatingPointType;

    SimpleFunctor(floatType cutoff)
        : autopas::Functor<Particle, ParticleCell>(cutoff), _cutoffSquare(cutoff * cutoff){};

    bool isRelevantForTuning() override { return true; }

    bool allowsNewton3() override { return true; }

    bool allowsNonNewton3() override { return false; }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      auto coordsI = i.getR();
      auto coordsJ = j.getR();

      std::array<double, 3> dr = autopas::ArrayMath::sub(coordsI, coordsJ);
      double dr2 = autopas::ArrayMath::dot(dr, dr);
      // in a grid with separation 1 this includes all neighbors with a Chebyshev distance of 1
      if (dr2 > _cutoffSquare) return;

      std::array<double, 3> f = {};

      for (int dim = 0; dim < 3; ++dim) {
        if (coordsI[dim] < coordsJ[dim]) {
          f[dim] = -1;
        } else if (coordsI[dim] > coordsJ[dim]) {
          f[dim] = 1;
        } else {
          f[dim] = 0;
        }
      }

      i.addF(f);
      j.subF(f);
    }

   private:
    const double _cutoffSquare;
  };
};
