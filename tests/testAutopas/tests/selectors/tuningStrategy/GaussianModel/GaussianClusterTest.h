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

namespace GaussianClusterTest {

class GaussianClusterTest : public AutoPasTestBase {
 protected:
  /**
   * Test optimization of a 1D discrete and 2D continuous function. The test
   * should find a value close to the optimum with a limited number of evaluation
   * of the black-box function.
   * @param functions vector of functions each corresponding to a discrete value
   * @param targetDiscrete optimal discrete value
   * @param targetContinuous optimal continuous value
   * @param precision Allowed difference between target and what is predicted.
   * @param domain domain of continuous part
   * @param acquisitionFunctionOption
   * @param visualize if true, the acquisition map is printed to std::cout
   */
  template <class NumberSetType>
  void test2DFunctions(const std::vector<std::function<double(double, double)>> &functions,
                       autopas::GaussianModelTypes::NeighbourFunction neighboursFun, int targetDiscrete,
                       const Eigen::VectorXd &targetContinuous, double precision,
                       const std::pair<NumberSetType, NumberSetType> &domain,
                       autopas::AcquisitionFunctionOption acquisitionFunctionOption, bool visualize) {
    autopas::Random rng(42);  // random generator

    constexpr size_t numEvidence = 10;     // number of samples allowed to make
    constexpr size_t lhsNumSamples = 850;  // number of samples to find max of acquisition function

    autopas::GaussianCluster gc({static_cast<int>(functions.size())}, 2,
                                autopas::GaussianCluster::WeightFunction::wasserstein2, 0.001, rng);

    size_t idEvidence = 0;

    auto evidenceContinuous = autopas::utils::Math::makeVectorXd({0, 0});
    for (idEvidence = 0; idEvidence < functions.size(); ++idEvidence) {
      auto evidenceDiscrete = autopas::utils::Math::makeVectorXi({static_cast<int>(idEvidence)});
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

      auto [amDiscrete, amContinuous] =
          gc.sampleAcquisitionMax(acquisitionFunctionOption, neighboursFun, lhsSamples).first;
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
                        autopas::GaussianModelTypes::NeighbourFunction neighboursFun,
                        autopas::AcquisitionFunctionOption af);

  /**
   * Print the input-output pair into std::cout
   * @param vecDiscrete
   * @param vecContinuous
   * @param out
   * @param evidenceNum
   */
  static void printEvidence(const Eigen::VectorXi &vecDiscrete, const Eigen::VectorXd &vecContinuous, double out,
                            size_t evidenceNum);
};

}  // end namespace GaussianClusterTest
