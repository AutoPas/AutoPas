/**
 * @file SimilarityFunctionsTest.cpp
 * @author S. Newcome
 * @date 06/06/2024
 */

#include <gtest/gtest.h>

#include <vector>

#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/utils/SimilarityFunctions.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Test calculateHomogeneityAndMaxDensity by
 * - creating a container with 88 particles, such that there are 8 bins
 * - using 88 particles, gives 11 particles per bin, testing handling of non-perfectly-divisible-by-10 numbers of
 * particles.
 * - the container has a cuboid shape to test handling of non cube shapes and is offset from the origin to test this
 *   scenario which is important e.g. if MPI is used.
 * - particles are placed evenly between bins
 * - tests density is correct and homogeneity is 0
 * - moves a particle between bins and tests again
 * - contiane
 */
TEST(SimilarityFunctionsTest, testCalculateHomogeneityAndMaxDensity) {
  const std::array<double, 3> domainMin{1, 1, 0};
  const std::array<double, 3> domainMax{9, 9, 4};

  // We arbitrarily use direct sum. Any container that implements ParticleContainerInterface should work.
  // cutoff = 3., skinPerTimestep = 0.1, verletRebuildFrequency = 10 are all abitrary
  autopas::DirectSum<Molecule> dsContainer{domainMin, domainMax, 3., 0.1, 10};

  const auto domainVol = (domainMax[0] - domainMin[0]) * (domainMax[1] - domainMin[1]) * (domainMax[2] - domainMin[2]);

  for (int i = 0; i < 88; ++i) {
    const size_t boxIndex1D = i / 11;
    const auto boxIndex3D = autopas::utils::ThreeDimensionalMapping::oneToThreeD(boxIndex1D, {2, 2, 2});
    const std::array<double, 3> position = {3. + (double)boxIndex3D[0] * 4., 3. + (double)boxIndex3D[1] * 4.,
                                            1. + (double)boxIndex3D[2] * 2.};

    ASSERT_LT(boxIndex3D[0], 2) << "Only two bins per direction expected\n";
    ASSERT_LT(boxIndex3D[1], 2) << "Only two bins per direction expected\n";
    ASSERT_LT(boxIndex3D[2], 2) << "Only two bins per direction expected\n";

    Molecule m(position, {0, 0, 0}, i);
    dsContainer.addParticle(m);
  }

  const double binVol = 4. * 4. * 2.;
  const double expectedMeanDensity = 88. / domainVol;
  const int numBins = 8;

  {
    const double expectedHomogeneity = 0.;
    const double expectedMaxDensity = expectedMeanDensity;

    const auto [actualHomogeneity, actualMaxDensity] = autopas::utils::calculateHomogeneityAndMaxDensity(dsContainer);

    EXPECT_NEAR(actualHomogeneity, expectedHomogeneity, 1e-12);
    EXPECT_NEAR(actualMaxDensity, expectedMaxDensity, 1e-12);
  }

  // Move particle 0

  for (auto iter = dsContainer.begin(); iter.isValid(); ++iter) {
    if (iter->getID() == 0) {
      iter->setR({7, 3, 1});  // Shifted to neighboring bin
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
    const double expectedMaxDensity = 12. / binVol;

    const auto [actualHomogeneity, actualMaxDensity] = autopas::utils::calculateHomogeneityAndMaxDensity(dsContainer);

    EXPECT_NEAR(actualHomogeneity, expectedHomogeneity, 1e-12);
    EXPECT_NEAR(actualMaxDensity, expectedMaxDensity, 1e-12);
  }
}