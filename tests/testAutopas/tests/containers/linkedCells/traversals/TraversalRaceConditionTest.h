/**
 * @file TraversalRaceConditionTest.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include <array>

#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

class TraversalRaceConditionTest : public AutoPasTestBase {
 public:
  TraversalRaceConditionTest() = default;

  ~TraversalRaceConditionTest() override = default;

  /*
   * Simple AoS only functor which repulses paritcles from each other with a
   * constant force of 1.
   */
  class SimpleFunctor : public autopas::Functor<Particle> {
   public:
    using SoAArraysType = Particle::SoAArraysType;
    using floatType = double;

    SimpleFunctor(floatType cutoff)
        : autopas::Functor<Particle>(cutoff), _cutoffSquare(cutoff * cutoff){};

    bool isRelevantForTuning() override { return true; }

    bool allowsNewton3() override { return true; }

    bool allowsNonNewton3() override { return false; }

    bool isAppropriateClusterSize(unsigned int clusterSize,
                                  autopas::DataLayoutOption::Value dataLayout) const override {
      return dataLayout == autopas::DataLayoutOption::aos;  // this functor supports clusters only for aos!
    }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      if (i.isDummy() or j.isDummy()) {
        return;
      }
      auto coordsI = i.getR();
      auto coordsJ = j.getR();

      std::array<double, 3> dr = autopas::utils::ArrayMath::sub(coordsI, coordsJ);
      double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
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
