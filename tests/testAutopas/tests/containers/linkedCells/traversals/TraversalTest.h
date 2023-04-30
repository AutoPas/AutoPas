/**
 * @file TraversalTest.h
 * @author C. Menges
 * @date 16.04.2019
 */

#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopasTools/generators/GridGenerator.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Test to check if all traversals consider all particles within cutoff
 */
class TraversalTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<autopas::TraversalOption, bool, autopas::LoadEstimatorOption>> {
 public:
  void SetUp() override {}

  void TearDown() override {}

  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      auto traversal(std::get<0>(inputTuple).to_string());
      auto loadEstimator(std::get<2>(inputTuple).to_string());
      // replace all '-' with '_', otherwise the test name is invalid
      std::replace(traversal.begin(), traversal.end(), '-', '_');
      std::replace(loadEstimator.begin(), loadEstimator.end(), '-', '_');
      return traversal + "_" + loadEstimator + "_" + (std::get<1>(inputTuple) ? "N3on" : "N3off");
    }
  };

  class CountFunctor : public autopas::Functor<Particle, CountFunctor> {
   public:
    using SoAArraysType = Particle::SoAArraysType;
    using ParticleCell = FPCell;
    using floatType = double;

    CountFunctor(floatType cutoff) : autopas::Functor<Particle, CountFunctor>(cutoff), _cutoffSquare(cutoff * cutoff){};

    bool isRelevantForTuning() override { return true; }

    bool allowsNewton3() override { return true; }

    bool allowsNonNewton3() override { return true; }

    bool allowsMixedNewton3() override { return false; }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      using namespace autopas::utils::ArrayMath::literals;

      if (i.isDummy() or j.isDummy()) {
        return;
      }
      const auto coordsI = i.getR();
      const auto coordsJ = j.getR();

      std::array<double, 3> dr = coordsI - coordsJ;
      const double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

      if (dr2 > _cutoffSquare) return;

      countFunc(i.getID());

      if (newton3) {
        countFunc(j.getID());
      }
    }

    // countFunc isn't a real function. It's just a declaration which is used to count the number of interactions.
    MOCK_METHOD1(countFunc, void(unsigned long id));

   private:
    floatType _cutoffSquare;
  };
};
