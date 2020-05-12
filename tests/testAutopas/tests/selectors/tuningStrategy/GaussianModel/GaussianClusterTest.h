/**
 * @file GaussianClusterTest.h
 * @author Jan Nguyen
 * @date 21.04.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "Eigen/Dense"
#include "autopas/selectors/FeatureVector.h"
#include "autopas/selectors/tuningStrategy/GaussianModel/GaussianCluster.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"

class GaussianClusterTest : public AutoPasTestBase {
 protected:
  /**
   *
   * @param functions vector of functions
   * @param targetDiscrete
   * @param targetContinuous
   * @param precision Allowed difference between target and what is predicted.
   * @param domain
   * @param acquisitionFunctionOption
   * @param visualize if true, the acquisition map is printed to std::cout
   */
  template <class NumberSetType>
  void test2DFunctions(const std::vector<std::function<double(double, double)>> &functions,
                       std::function<std::vector<Eigen::VectorXi>(Eigen::VectorXi)> neighboursFun, int targetDiscrete,
                       const Eigen::VectorXd &targetContinuous, double precision,
                       const std::pair<NumberSetType, NumberSetType> &domain,
                       autopas::AcquisitionFunctionOption acquisitionFunctionOption, bool visualize) {
    autopas::Random rng(42);  // random generator

    constexpr size_t numEvidence = 10;     // number of samples allowed to make
    constexpr size_t lhsNumSamples = 850;  // number of samples to find max of acquisition function

    autopas::GaussianCluster gc({static_cast<int>(functions.size())}, 2, 0.001, rng);

    size_t idEvidence = 0;

    Eigen::VectorXd evidenceContinuous(2);
    evidenceContinuous << 0, 0;
    for (idEvidence = 0; idEvidence < functions.size(); ++idEvidence) {
      Eigen::VectorXi evidenceDiscrete(1);
      evidenceDiscrete << idEvidence;

      auto evidenceOut = functions[evidenceDiscrete[0]](evidenceContinuous[0], evidenceContinuous[1]);

      if (visualize) {
        printEvidence(evidenceDiscrete, evidenceContinuous, evidenceOut, idEvidence);
      }
      gc.addEvidence(evidenceDiscrete, evidenceContinuous, evidenceOut);
    }

    for (; idEvidence < numEvidence; ++idEvidence) {
      // create lhs samples
      std::vector<Eigen::VectorXd> lhsSamples;
      lhsSamples.reserve(lhsNumSamples);

      auto xSamples = domain.first.uniformSample(lhsNumSamples, rng);
      auto ySamples = domain.second.uniformSample(lhsNumSamples, rng);
      for (size_t idSample = 0; idSample < lhsNumSamples; ++idSample) {
        Eigen::Vector2d sample(xSamples[idSample], ySamples[idSample]);
        lhsSamples.emplace_back(sample);
      }

      // calculate all acquisitions
      auto acquisitions = gc.sampleAcquisition(acquisitionFunctionOption, neighboursFun, lhsSamples);

      const auto &[amDiscrete, amContinuous] = acquisitions.back();
      double amOut = functions[amDiscrete[0]](amContinuous[0], amContinuous[1]);

      if (visualize) {
        printMaps(20, 20, domain.first, domain.second, gc, neighboursFun, acquisitionFunctionOption);
        printEvidence(amDiscrete, amContinuous, amOut, idEvidence);
      }

      gc.addEvidence(amDiscrete, amContinuous, amOut);
    }

    auto [predictedMaxDiscrete, predictedMaxContinuous] = gc.getEvidenceMax();

    // check if prediction is near target
    double predMax = functions[predictedMaxDiscrete[0]](predictedMaxContinuous[0], predictedMaxContinuous[1]);
    double realMax = functions[targetDiscrete](targetContinuous[0], targetContinuous[1]);
    EXPECT_NEAR(predMax, realMax, precision);
  }

  /**
   * Print a xChunks x yChunks aquisition map for each functor.
   * @param xChunks width of maps
   * @param yChunks height of maps
   * @param domainX x domain of GaussianProcess
   * @param domainY y domain of GaussianProcess
   * @param gc GaussianCluster to calculate acquisiton
   * @param af acquisition function
   */
  static void printMaps(int xChunks, int yChunks, const autopas::NumberSet<double> &domainX,
                        const autopas::NumberSet<double> &domainY, const autopas::GaussianCluster &gc,
                        std::function<std::vector<Eigen::VectorXi>(Eigen::VectorXi)> neighboursFun,
                        autopas::AcquisitionFunctionOption af);

  /**
   * Print the input-output pair into std::cout
   * @param vecDiscrete
   * @param vecContinuous
   * @param out
   * @param evidenceNum
   */
  static void printEvidence(Eigen::VectorXi vecDiscrete, Eigen::VectorXd vecContinuous, double out, size_t evidenceNum);
};
