/**
* @file LiveInfoTest.cpp
* @author S. Newcome
* @date 02/04/2025
*/

#include "LiveInfoTest.h"
#include "autopas/tuning/tuningStrategy/LiveInfo.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ParticleBinStructure.h"
#include "autopas/utils/ArrayMath.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Tests Live Info in the following ways:
 * - Adds 50 owned, halo, and dummy particles to a 10x10x10 domain.
 * - Checks that the statistics calculated within LiveInfo directly are correct compared to expected "hand-calculated"
 *   values.
 * - Checks that the LiveInfo statistics calculated through ParticleBinStructure within LiveInfo match directly using
 *   ParticleBinStructure (that these statistics are correct is not the responsibility of this test).
 */
TEST_F(LiveInfoTest, CorrectInfoMappings) {
  using namespace autopas::utils::ArrayMath;
  using namespace literals;


  autopas::LiveInfo liveInfo;

  constexpr double cutoff{2.0};
  constexpr double skin{0.1};
  constexpr std::array<double, 3> boxMin{0., 0., 0.};
  constexpr std::array<double, 3> boxMax{10., 10., 10.};
  constexpr auto boxSize = boxMax - boxMin;
  constexpr size_t rebuildFrequency{10};
  constexpr size_t numParticlesOwned{50};
  constexpr size_t numParticlesHalo{30};
  constexpr size_t numParticlesDummy{10};

  auto container = autopas::AutoPas<ParticleFP64>();

  container.setCutoff(cutoff);
  container.setVerletSkin(skin);
  container.setBoxMin(boxMin);
  container.setBoxMax(boxMax);
  container.setVerletRebuildFrequency(rebuildFrequency);

  container.init();


  // Create comparison bins
  const auto numCellBinsPerDim =  castedFloor<size_t>(boxSize / (cutoff + skin));
  const auto numCells = numCellBinsPerDim[0] * numCellBinsPerDim[1] * numCellBinsPerDim[2];
  const auto cellBinDims = boxSize / staticCastArray<double>(numCellBinsPerDim);
  auto cellBins = autopas::utils::ParticleBinStructure(numCellBinsPerDim, cellBinDims, boxMin, boxMax, cutoff);

  constexpr auto blurredBinDims = boxSize / std::array<double, 3>{3., 3., 3.};
  auto blurredBins = autopas::utils::ParticleBinStructure({3, 3, 3}, blurredBinDims, boxMin, boxMax, cutoff);


  // Add particles to the container and bin them
  size_t id{0};
  for (size_t i = 0; i < numParticlesOwned; ++i, ++id) {
    ParticleFP64 p({5., 5., 5.}, {0., 0., 0.}, id);
    container.addParticle(p);
    cellBins.countParticle({5., 5., 5.});
    blurredBins.countParticle({5., 5., 5.});
  }
  for (size_t i = 0; i < numParticlesHalo; ++i, ++id) {
    ParticleFP64 p({-1., -1., -1.}, {0., 0., 0.}, id);
    container.addHaloParticle(p);
  }
  // Cannot directly add dummy particles, so will add the particles as owned, then find them and turn them into dummys.
  for (size_t i = 0; i < numParticlesDummy; ++i, ++id) {
    ParticleFP64 p({1., 1., 1.}, {0., 0., 0.}, id);
    container.addParticle(p);
  }
  for (auto p = container.begin(); p.isValid(); ++p) {
    if (p->getID() >= numParticlesOwned + numParticlesHalo) {
      p->setOwnershipState(autopas::OwnershipState::dummy);
    }
  }


  // Gather information from the container
  liveInfo.gather(container.begin(), container.getVerletRebuildFrequency(), container.getNumberOfParticles(), container.getBoxMin(),
    container.getBoxMax(), container.getCutoff(), container.getVerletSkin());

  cellBins.calculateStatistics();
  blurredBins.calculateStatistics();

  // Compare directly calculated LiveInfo statistics
  EXPECT_EQ(liveInfo.get<size_t>("numOwnedParticles"), numParticlesOwned);
  EXPECT_EQ(liveInfo.get<size_t>("numHaloParticles"), numParticlesHalo);
  EXPECT_EQ(liveInfo.get<double>("cutoff"), cutoff);
  EXPECT_EQ(liveInfo.get<double>("skin"), skin);
  EXPECT_EQ(liveInfo.get<size_t>("rebuildFrequency"), rebuildFrequency);
  EXPECT_EQ(liveInfo.get<size_t>("particleSize"), sizeof(ParticleFP64));
  EXPECT_EQ(liveInfo.get<size_t>("threadCount"), static_cast<size_t>(autopas::autopas_get_max_threads()));
  EXPECT_EQ(liveInfo.get<double>("domainSizeX"), boxSize[0]);
  EXPECT_EQ(liveInfo.get<double>("domainSizeY"), boxSize[1]);
  EXPECT_EQ(liveInfo.get<double>("domainSizeZ"), boxSize[2]);
  EXPECT_EQ(liveInfo.get<size_t>("numCells"), numCells);

  // Compare LiveInfo statistics calculated through ParticleBinStructure
  EXPECT_EQ(liveInfo.get<size_t>("numEmptyCells"), cellBins.getNumEmptyBins());
  EXPECT_EQ(liveInfo.get<size_t>("maxParticlesPerCell"), cellBins.getMaxParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("minParticlesPerCell"), cellBins.getMinParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("medianParticlesPerCell"), cellBins.getMedianParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("lowerQuartileParticlesPerCell"), cellBins.getLowerQuartileParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("upperQuartileParticlesPerCell"), cellBins.getUpperQuartileParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("meanParticlesPerCell"), cellBins.getMeanParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("relativeParticlesPerCellStdDev"), cellBins.getRelStdDevParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("particlesPerCellStdDev"), cellBins.getStdDevParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("estimatedNumNeighborInteractions"), cellBins.getEstimatedNumberOfNeighborInteractions());

  EXPECT_EQ(liveInfo.get<size_t>("maxParticlesPerBlurredBin"), blurredBins.getMaxParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("minParticlesPerBlurredBin"), blurredBins.getMinParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("medianParticlesPerBlurredBin"), blurredBins.getMedianParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("lowerQuartileParticlesPerBlurredBin"), blurredBins.getLowerQuartileParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("upperQuartileParticlesPerBlurredBin"), blurredBins.getUpperQuartileParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("meanParticlesPerBlurredBin"), blurredBins.getMeanParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("relativeParticlesPerBlurredBinStdDev"), blurredBins.getRelStdDevParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("particlesPerBlurredBinStdDev"), blurredBins.getStdDevParticlesPerBin());

}

/**
 * Tests that no LiveInfo Statistics are missing from the above test.
 *
 * There should be
 * - 11 directly calculated statistics
 * - 10 cellBin statistics
 * - 8 blurredBin statistics
 *
 * = 29
 */
TEST_F(LiveInfoTest, NoMissingLiveInfoStatisticsUnitTests) {
  autopas::LiveInfo liveInfo;
  auto container = autopas::AutoPas<ParticleFP64>();

  container.setCutoff(2.);
  container.setVerletSkin(0.1);
  container.setBoxMin({0., 0., 0.});
  container.setBoxMax({10., 10., 10.});
  container.setVerletRebuildFrequency(5);

  container.init();

  ParticleFP64 p({5., 5., 5.}, {0., 0., 0.}, 0);

  container.addParticle(p);

  // Gather information from the container
  liveInfo.gather(container.begin(), container.getVerletRebuildFrequency(), container.getNumberOfParticles(), container.getBoxMin(),
    container.getBoxMax(), container.getCutoff(), container.getVerletSkin());

  EXPECT_EQ(size(liveInfo.get()), 29) << "LiveInfoTest is missing a test for one of the LiveInfo statistics, or"
                                         " includes a test for a statistic that no longer exists.";
}