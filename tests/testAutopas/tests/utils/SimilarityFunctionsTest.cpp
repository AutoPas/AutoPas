/**
* @file SimilarityFunctionsTest.cpp
* @author S. Newcome
* @date 06/06/2024
 */

#include <algorithm>
#include <iterator>
#include <vector>

#include <gtest/gtest.h>
#include "autopas/AutoPasDecl.h"
#include "testingHelpers/commonTypedefs.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/SimilarityFunctions.h"

/**
 * Test calculateHomogeneityAndMaxDensity by
 * - creating a container with 80 particles, such that there are 8 bins
 * - the container has a cuboid shape to test handling of non cube shapes and is offset from the origin to test this
 *   scenario which is important e.g. if MPI is used.
 * - particles are placed evenly between bins
 * - tests density is correct and homogeneity is 0
 * - moves a particle between bins and tests again
 * - contiane
 */
TEST(SimilarityFunctionsTest, testCalculateHomogeneityAndMaxDensity) {
  autopas::AutoPas<Molecule> autoPas;

  const std::array<double, 3> domainMin{1,1,0};
  const std::array<double, 3> domainMax{9,9,8};
  const auto domainVol = (domainMax[0] - domainMin[0]) * (domainMax[1] - domainMin[1]) * (domainMax[2] - domainMin[2]);

  autoPas.setBoxMin(domainMin);
  autoPas.setBoxMax(domainMax);

  autoPas.init();

  for (int i = 0; i < 80; ++i) {
    const auto boxIndex1D = i / 10;
    const auto boxIndex3D = autopas::utils::ThreeDimensionalMapping::oneToThreeD(boxIndex1D, {2,2,2});
    const std::array<double, 3> position = {3. + (double)boxIndex3D[0] * 4., 3. + (double)boxIndex3D[1] * 4., 2. + (double)boxIndex3D[2] * 4.};

    Molecule m(position, {0,0,0}, i);
    autoPas.addParticle(m);
  }

  const double binVol = 4. * 4. * 4.;
  const double expectedMeanDensity = 80. / domainVol;
  const int numBins = 8;

  {
    const double expectedHomogeneity = 0.;
    const double expectedMaxDensity = expectedMeanDensity;

    const auto [actualHomogeneity, actualMaxDensity] = autopas::utils::calculateHomogeneityAndMaxDensity(autoPas, autoPas.getBoxMin(), autoPas.getBoxMax());

    EXPECT_NEAR(actualHomogeneity, expectedHomogeneity, 1e-12);
    EXPECT_NEAR(actualMaxDensity, expectedMaxDensity, 1e-12);
  }

  // Move particle 0

  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 0) {
      iter->setR({7,3,2}); // Shifted to neighboring bin
      break;
    }
  }

  {
    // Most bins are still at the mean density but two have different densities
    const double densityDiff = 1. / binVol;

    // Both denser and less dense bins have the same varaince contribution
    const auto varianceContributionOfDiffBins = densityDiff * densityDiff / (double)numBins;

    // Add contributions from both bins.
    const auto densityVariance = 2. * varianceContributionOfDiffBins;

    const double expectedHomogeneity = std::sqrt(densityVariance);
    const double expectedMaxDensity = 11. / binVol;

    const auto [actualHomogeneity, actualMaxDensity] = autopas::utils::calculateHomogeneityAndMaxDensity(autoPas, autoPas.getBoxMin(), autoPas.getBoxMax());

    EXPECT_NEAR(actualHomogeneity, expectedHomogeneity, 1e-12);
    EXPECT_NEAR(actualMaxDensity, expectedMaxDensity, 1e-12);
  }

}