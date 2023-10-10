/**
 * @file GaussianProcessTest.h
 * @author Jan Nguyen
 * @date 12.06.19
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "Eigen/Dense"
#include "autopas/tuning/tuningStrategy/GaussianModel/GaussianProcess.h"
#include "autopas/tuning/utils/FeatureVector.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"

class GaussianProcessTest : public AutoPasTestBase {
 protected:
  /**
   *
   * @param function
   * @param target
   * @param precision Allowed difference between target and what is predicted.
   * @param domain
   * @param acquisitionFunctionOption
   * @param visualize if true, the acquisition map is printed to std::cout
   */
  template <class NumberSetType>
  void test2DFunction(const std::function<double(double, double)> &function, const Eigen::VectorXd &target,
                      double precision, const std::pair<NumberSetType, NumberSetType> &domain,
                      autopas::AcquisitionFunctionOption acquisitionFunctionOption, bool visualize) {
    autopas::Random rng(32);  // random generator

    constexpr size_t numEvidence = 8;      // number of samples allowed to make
    constexpr size_t lhsNumSamples = 900;  // number of samples to find max of acquisition function

    autopas::GaussianProcess gp(2, 0.001, rng);

    // add first evidence
    auto firstEvidence = autopas::utils::Math::makeVectorXd({0, 0});
    gp.addEvidence(firstEvidence, function(0, 0), true);

    for (unsigned idEvidence = 1; idEvidence < numEvidence; ++idEvidence) {
      // create lhs samples
      std::vector<Eigen::VectorXd> lhsSamples;
      lhsSamples.reserve(lhsNumSamples);

      auto xSamples = domain.first.uniformSample(lhsNumSamples, rng);
      auto ySamples = domain.second.uniformSample(lhsNumSamples, rng);
      for (size_t idSample = 0; idSample < lhsNumSamples; ++idSample) {
        Eigen::Vector2d sample(xSamples[idSample], ySamples[idSample]);
        lhsSamples.emplace_back(sample);
      }

      // sample acquisition function
      auto am = gp.sampleAquisitionMax(acquisitionFunctionOption, lhsSamples);
      double amOut = function(am[0], am[1]);

      if (visualize) {
        printMap(20, 20, domain.first, domain.second, gp, acquisitionFunctionOption, 1.);
        std::cout << "Acq max: " << std::endl << am << std::endl;
        std::cout << "Got: " << amOut << std::endl;
      }
      gp.addEvidence(am, amOut, true);
    }

    auto predictedTarget = gp.getEvidenceMax();

    // check if prediction is near target
    double predMax = function(predictedTarget[0], predictedTarget[1]);
    double realMax = function(target[0], target[1]);
    EXPECT_NEAR(predMax, realMax, precision);
  }

  /**
   * Print a xChunks x yChunks mean and aquisition map.
   * @param xChunks width of maps
   * @param yChunks height of maps
   * @param domainX x domain of GaussianProcess
   * @param domainY y domain of GaussianProcess
   * @param gp GaussianProcess to calculate mean and acquisiton
   * @param af acquisition function
   * @param colorFactor Multiplies the mean by this factor to gain the color
   */
  static void printMap(int xChunks, int yChunks, const autopas::NumberSet<double> &domainX,
                       const autopas::NumberSet<double> &domainY, const autopas::GaussianProcess &gp,
                       autopas::AcquisitionFunctionOption af, double colorFactor);
};
