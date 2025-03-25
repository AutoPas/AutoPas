/**
 * @file ParticleBinStructureTests.cpp
 * @author S. Newcome
 * @date 24/03/2025
 */

#include "autopas/utils/ParticleBinStructure.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

#include <gtest/gtest.h>
#include <random>


/**
 * Adds six particles to a bin structure (or attempts to) and checks that the total and individual particle counts are
 * correct.
 */
TEST(ParticleBinStructureTest, testCountParticle) {
  using autopas::utils::ArrayMath::isEqual;

  // Construct a bin structure with 2 bins in each dimension, each bin is 1x1x1
  auto particleBinStructure = autopas::utils::ParticleBinStructure({2, 2, 2}, {1., 1., 1.}, {0., 0., 0.}, {2. ,2. ,2.}, 1.);

  // Check counts are 0
  EXPECT_EQ(particleBinStructure.getTotalParticleCount(), 0);
  EXPECT_EQ(particleBinStructure.getParticleCounts(), std::vector<size_t>({0, 0, 0, 0, 0, 0, 0, 0}));

  // Add a particle to the first bin
  particleBinStructure.countParticle({0.5, 0.5, 0.5});
  EXPECT_EQ(particleBinStructure.getTotalParticleCount(), 1);
  EXPECT_EQ(particleBinStructure.getParticleCounts(), std::vector<size_t>({1, 0, 0, 0, 0, 0, 0, 0}));

  // Add a particle to the second bin
  particleBinStructure.countParticle({1.5, 0.5, 0.5});
  EXPECT_EQ(particleBinStructure.getTotalParticleCount(), 2);
  EXPECT_EQ(particleBinStructure.getParticleCounts(), std::vector<size_t>({1, 1, 0, 0, 0, 0, 0, 0}));

  // Add a particle to the first bin
  particleBinStructure.countParticle({0.5, 0.5, 0.5});
  EXPECT_EQ(particleBinStructure.getTotalParticleCount(), 3);
  EXPECT_EQ(particleBinStructure.getParticleCounts(), std::vector<size_t>({2, 1, 0, 0, 0, 0, 0, 0}));

  // Add a particle on the edge of the structure to check that it is clamped to the last bin
  particleBinStructure.countParticle({2., 2., 2.});
  EXPECT_EQ(particleBinStructure.getTotalParticleCount(), 4);
  EXPECT_EQ(particleBinStructure.getParticleCounts(), std::vector<size_t>({2, 1, 0, 0, 0, 0, 0, 1}));

  // Add a particle outside the structure to check that it is ignored
  particleBinStructure.countParticle({3., 3., 3.});
  EXPECT_EQ(particleBinStructure.getTotalParticleCount(), 4);
  EXPECT_EQ(particleBinStructure.getParticleCounts(), std::vector<size_t>({2, 1, 0, 0, 0, 0, 0, 1}));

  // Add a particle on an internal bin boundary. We don't care which bin it goes into, just that it goes into one.
  particleBinStructure.countParticle({1., 0.5, 0.5});
  EXPECT_EQ(particleBinStructure.getTotalParticleCount(), 5);
  EXPECT_TRUE(isEqual(particleBinStructure.getParticleCounts(), std::vector<size_t>({3, 1, 0, 0, 0, 0, 0, 1})) or
              isEqual(particleBinStructure.getParticleCounts(),std::vector<size_t>({2, 2, 0, 0, 0, 0, 0, 1})));

}

/**
 * Test for calculateStatistics. This is a basic test intended to check that the function runs and to check basic
 * behavior without comparing to reference values.
 */
TEST(ParticleBinStructureTest, testCalculateStatisticsBasic) {
  using namespace autopas::utils;
  using namespace ArrayMath::literals;

  const std::array<size_t, 3> numBinsPerDim = {5, 5, 5};
  const size_t numBins = numBinsPerDim[0] * numBinsPerDim[1] * numBinsPerDim[2];
  const std::array<double, 3> binLength = {1.2, 1., 1.};
  const std::array<double, 3> boxMax = ArrayMath::staticCastArray<double>(numBinsPerDim) * binLength;
  auto particleBinStructure = autopas::utils::ParticleBinStructure(numBinsPerDim, binLength, {0., 0., 0.}, boxMax, 1.);

  // Calculate reference estimatedHitRate
  const auto volCutoffSphere = 4. / 3. * M_PI * 1. * 1. * 1.;
  const auto volBin = 1.2 * 1. * 1.;
  const auto potentialInteractionVol = volBin * 27.;
  const auto estimatedHitRate = volCutoffSphere / potentialInteractionVol;

  particleBinStructure.calculateStatistics();

  // Check that all statistics correctly reflect there being no particles in the structure
  EXPECT_DOUBLE_EQ(particleBinStructure.getMeanParticlesPerBin(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getStdDevParticlesPerBin(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getRelStdDevParticlesPerBin(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getMeanDensity(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getStdDevDensity(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getMaxDensity(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getEstimatedNumberOfNeighborInteractions(), 0.);
  EXPECT_EQ(particleBinStructure.getMinParticlesPerBin(), 0);
  EXPECT_EQ(particleBinStructure.getMaxParticlesPerBin(), 0);
  EXPECT_EQ(particleBinStructure.getLowerQuartileParticlesPerBin(), 0);
  EXPECT_EQ(particleBinStructure.getUpperQuartileParticlesPerBin(), 0);
  EXPECT_EQ(particleBinStructure.getMedianParticlesPerBin(), 0);
  EXPECT_EQ(particleBinStructure.getNumEmptyBins(), numBins);

  // Add 10 particles to each bin
  for (size_t i = 0; i < numBins; ++i) {
    const auto binIndex3D = ThreeDimensionalMapping::oneToThreeD(i, numBinsPerDim);
    const auto position = ArrayMath::staticCastArray<double>(binIndex3D) * binLength + binLength / 2.;
    for (size_t j = 0; j < 10; ++j) {
      particleBinStructure.countParticle(position);
    }
  }
  particleBinStructure.calculateStatistics();

  // Check that all statistics correctly reflect an equal number of particles in each bin (and that this number is 10).
  EXPECT_DOUBLE_EQ(particleBinStructure.getMeanParticlesPerBin(), 10.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getStdDevParticlesPerBin(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getRelStdDevParticlesPerBin(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getMeanDensity(), 10. / (binLength[0] * binLength[1] * binLength[2]));
  EXPECT_DOUBLE_EQ(particleBinStructure.getStdDevDensity(), 0.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getMaxDensity(), 10.);
  EXPECT_DOUBLE_EQ(particleBinStructure.getEstimatedNumberOfNeighborInteractions(), 10. * (10. * 27.) * estimatedHitRate - 10.);
  EXPECT_EQ(particleBinStructure.getMinParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getMaxParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getLowerQuartileParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getUpperQuartileParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getMedianParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getNumEmptyBins(), 0);

  // Add another particle. Check that all statistics increase (specifically that deviation based statistics are non-zero.)
  // and where reference values are available without calculation, check that they are correct.
  particleBinStructure.countParticle({5., 5., 5.});
  particleBinStructure.calculateStatistics();

  EXPECT_DOUBLE_EQ(particleBinStructure.getMeanParticlesPerBin(), 10.1);
  EXPECT_GT(particleBinStructure.getStdDevParticlesPerBin(), 0.);
  EXPECT_GT(particleBinStructure.getRelStdDevParticlesPerBin(), 0.);
  EXPECT_GT(particleBinStructure.getStdDevDensity(), 0.);
  EXPECT_GT(particleBinStructure.getMaxDensity(), 10.);
  EXPECT_GT(particleBinStructure.getEstimatedNumberOfNeighborInteractions(), 10. * (10. * 27.) * estimatedHitRate - 10.);
  EXPECT_EQ(particleBinStructure.getMinParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getMaxParticlesPerBin(), 11);
  EXPECT_EQ(particleBinStructure.getLowerQuartileParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getUpperQuartileParticlesPerBin(), 11);
  EXPECT_EQ(particleBinStructure.getMedianParticlesPerBin(), 10);
  EXPECT_EQ(particleBinStructure.getNumEmptyBins(), 0);
}

/**
 * A test to check that the statistics calculated by ParticleBinStructure are the same as those calculated by a
 * (deterministic) random reference.
 */
TEST(ParticleBinStructureTest, calculateStatisticsVsReference) {
  using namespace autopas::utils;
  using namespace ArrayMath::literals;

  const std::array<size_t, 3> numBinsPerDim = {5, 5, 5};
  const size_t numBins = numBinsPerDim[0] * numBinsPerDim[1] * numBinsPerDim[2];
  const std::array<double, 3> binLength = {1.2, 1., 1.};
  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = ArrayMath::staticCastArray<double>(numBinsPerDim) * binLength;
  auto particleBinStructure = autopas::utils::ParticleBinStructure(numBinsPerDim, binLength, boxMin, boxMax, 1.);

  const auto volBin = binLength[0] * binLength[1] * binLength[2];


  // Generate random distribution and count particles

  // Seed for deterministic random number generation
  std::mt19937 gen(42);
  std::uniform_real_distribution<> distributionX(boxMin[0], boxMax[0]);
  std::uniform_real_distribution<> distributionY(boxMin[1], boxMax[1]);
  std::uniform_real_distribution<> distributionZ(boxMin[2], boxMax[2]);

  // Fill the structure with random particle positions
  const size_t numParticles = 1000;
  for (size_t i = 0; i < numParticles; ++i) {
    std::array<double, 3> particlePosition = {distributionX(gen), distributionY(gen), distributionZ(gen)};
    particleBinStructure.countParticle(particlePosition);
  }

  // Calculate statistics
  particleBinStructure.calculateStatistics();

  // Calculate reference statistics and compare
  const auto particleCounts = particleBinStructure.getParticleCounts();

  // Mean
  const auto expectedMean = static_cast<double>(numParticles) / static_cast<double>(numBins);
  EXPECT_DOUBLE_EQ(particleBinStructure.getMeanParticlesPerBin(), expectedMean);

  // Std. Dev.
  const auto expectedStdDev = [&]() {
    double sum = 0.;
    for (size_t i = 0; i < numBins; ++i) {
      const auto diff = static_cast<double>(particleCounts[i]) - expectedMean;
      sum += diff * diff;
    }
    return std::sqrt(sum / static_cast<double>(numBins));
  }();
  EXPECT_DOUBLE_EQ(particleBinStructure.getStdDevParticlesPerBin(), expectedStdDev);

  // Rel. Std. Dev.
  const auto expectedRelStdDev = expectedStdDev / expectedMean;
  EXPECT_DOUBLE_EQ(particleBinStructure.getRelStdDevParticlesPerBin(), expectedRelStdDev);

  // Max
  const auto expectedMax = *std::max_element(particleCounts.begin(), particleCounts.end());
  EXPECT_EQ(particleBinStructure.getMaxParticlesPerBin(), expectedMax);

  // Min
  const auto expectedMin = *std::min_element(particleCounts.begin(), particleCounts.end());
  EXPECT_EQ(particleBinStructure.getMinParticlesPerBin(), expectedMin);

  // Median
  const auto expectedMedian = [&]() {
    std::vector<size_t> sortedCounts(particleCounts);
    std::sort(sortedCounts.begin(), sortedCounts.end());
    return sortedCounts[sortedCounts.size() / 2];
  }();
  EXPECT_EQ(particleBinStructure.getMedianParticlesPerBin(), expectedMedian);

  // Lower Quartile
  const auto expectedLowerQuartile = [&]() {
    std::vector<size_t> sortedCounts(particleCounts);
    std::sort(sortedCounts.begin(), sortedCounts.end());
    return sortedCounts[sortedCounts.size() / 4];
  }();
  EXPECT_EQ(particleBinStructure.getLowerQuartileParticlesPerBin(), expectedLowerQuartile);

  // Upper Quartile
  const auto expectedUpperQuartile = [&]() {
    std::vector<size_t> sortedCounts(particleCounts);
    std::sort(sortedCounts.begin(), sortedCounts.end());
    return sortedCounts[3 * sortedCounts.size() / 4];
  }();
  EXPECT_EQ(particleBinStructure.getUpperQuartileParticlesPerBin(), expectedUpperQuartile);

  // Mean Density
  const auto expectedMeanDensity = static_cast<double>(numParticles) / (boxMax[0] * boxMax[1] * boxMax[2]);
  EXPECT_DOUBLE_EQ(particleBinStructure.getMeanDensity(), expectedMeanDensity);

  // Std. Dev. Density
  const auto expectedStdDevDensity = [&]() {
    double sum = 0.;
    for (size_t i = 0; i < numBins; ++i) {
      const auto diff = (static_cast<double>(particleCounts[i]) / volBin) - expectedMeanDensity;
      sum += diff * diff;
    }
    return std::sqrt(sum / static_cast<double>(numBins));
  }();
  EXPECT_DOUBLE_EQ(particleBinStructure.getStdDevDensity(), expectedStdDevDensity);

  // Max Density
  const auto expectedMaxDensity = static_cast<double>(expectedMax) / volBin;
  EXPECT_DOUBLE_EQ(particleBinStructure.getMaxDensity(), expectedMaxDensity);

  // Num Empty Bins
  const auto expectedNumEmptyBins = std::count(particleCounts.begin(), particleCounts.end(), 0);
  EXPECT_EQ(particleBinStructure.getNumEmptyBins(), expectedNumEmptyBins);

  // Estimated Number of Neighbor Interactions
  const auto volCutoffSphere = 4. / 3. * M_PI * 1. * 1. * 1.;
  const auto potentialInteractionVol = volBin * 27.;
  const auto estimatedHitRate = volCutoffSphere / potentialInteractionVol;
  const auto expectedEstimatedNumberOfNeighborInteractions = [&]() {
    double sum = 0.;
    for (size_t i = 0; i < numBins; ++i) {
      const auto particleCount = particleCounts[i];
      if (particleCount == 0) {
        continue;
      }
      sum += std::max(static_cast<double>(particleCount * (particleCount * 27)) * estimatedHitRate - particleCount, 0.);
    }
    return sum;
  }();
  EXPECT_DOUBLE_EQ(particleBinStructure.getEstimatedNumberOfNeighborInteractions(), expectedEstimatedNumberOfNeighborInteractions);
}

/**
 * Tests estimatedNumberOfNeighborInteractions for a single bin against the actual number given the assumptions hold.
 * Uses a random uniform distribution of particles with sufficiently enough particles to avoid disrepencies from particle placement.
 */
TEST(ParticleBinStructureTest, estimatedNumNeighIntVsActualWithAssumptions) {
  using namespace autopas::utils;
  using namespace ArrayMath::literals;

  const std::array<size_t, 3> numBinsPerDim = {1, 1, 1};
  const size_t numBins = numBinsPerDim[0] * numBinsPerDim[1] * numBinsPerDim[2];
  const std::array<double, 3> binLength = {1., 1., 1.};
  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = ArrayMath::staticCastArray<double>(numBinsPerDim) * binLength;
  auto particleBinStructure = autopas::utils::ParticleBinStructure(numBinsPerDim, binLength, boxMin, boxMax, 1.);

  const auto volBin = binLength[0] * binLength[1] * binLength[2];


  // Generate random distribution and count particles
  // Generated particles stretch beyond the domain (e.g. like a halo layer) to avoid discrepancies.

  // Seed for deterministic random number generation
  std::mt19937 gen(42);
  std::uniform_real_distribution<> distributionX(boxMin[0]-1., boxMax[0]+1.);
  std::uniform_real_distribution<> distributionY(boxMin[1]-1., boxMax[1]+1.);
  std::uniform_real_distribution<> distributionZ(boxMin[2]-1., boxMax[2]+1.);

  // Fill the structure with random particle positions
  const size_t numParticles = 5000;
  std::vector<std::array<double, 3>> particlePositions(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    std::array<double, 3> particlePosition = {distributionX(gen), distributionY(gen), distributionZ(gen)};
    particlePositions.push_back(particlePosition);
    particleBinStructure.countParticle(particlePosition);
  }

  // Calculate statistics
  particleBinStructure.calculateStatistics();

  size_t numNeighborInteractions = 0;
  for (size_t i = 0; i < numParticles; ++i) {
    for (size_t j = i + 1; j < numParticles; ++j) {
      // Ignore particle pairs that are both "halo" particles
      if ((particlePositions[i][0] < boxMin[0] or particlePositions[i][0] > boxMax[0] or
          particlePositions[i][1] < boxMin[1] or particlePositions[i][1] > boxMax[1] or
          particlePositions[i][2] < boxMin[2] or particlePositions[i][2] > boxMax[2]) and
          (particlePositions[j][0] < boxMin[0] or particlePositions[j][0] > boxMax[0] or
          particlePositions[j][1] < boxMin[1] or particlePositions[j][1] > boxMax[1] or
          particlePositions[j][2] < boxMin[2] or particlePositions[j][2] > boxMax[2])) {
        continue;
      }
      const auto displacement = particlePositions[i] - particlePositions[j];
      const auto distanceSquared = ArrayMath::dot(displacement, displacement);
      if (distanceSquared < 1.) {
        numNeighborInteractions++;
      }
    }
  }

  EXPECT_NEAR(particleBinStructure.getEstimatedNumberOfNeighborInteractions(), numNeighborInteractions, 10);
}