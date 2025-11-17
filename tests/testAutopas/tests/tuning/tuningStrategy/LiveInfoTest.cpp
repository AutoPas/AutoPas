/**
 * @file LiveInfoTest.cpp
 * @author S. Newcome
 * @date 02/04/2025
 */

#include "LiveInfoTest.h"

#include "autopas/AutoPas.h"
#include "autopas/tuning/tuningStrategy/LiveInfo.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ParticleBinStructure.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Tests that the parametrized constructor creates the correct info mappings.
 */
TEST_F(LiveInfoTest, ParametrizedConstructorCreatesCorrectInfo) {
  const std::vector<std::tuple<size_t, size_t, double, double, double, double, double, size_t, size_t, size_t, size_t,
                               size_t, size_t, size_t, size_t, size_t, size_t, double, double, double, size_t, double,
                               double, size_t, size_t, size_t, size_t, size_t, double, double, double>>
      parameters = {
          {10, 5,   2.0, 0.1,  10.0, 10.0, 10.0, 64, 4, 20, 8, 2, 1,   3,   2,   1,
           3,  2.0, 0.5, 0.25, 100,  0.8,  0.05, 4,  1, 2,  1, 3, 2.0, 0.3, 0.15},
          {50, 25,  3.0, 0.2, 20.0, 20.0, 20.0, 128, 8, 30, 16, 4, 2,   6,   4,  2,
           6,  3.0, 1.0, 0.5, 200,  1.6,  0.1,  8,   2, 4,  2,  6, 3.0, 0.6, 0.3},
          {25, 8,   0.750, 0.15, 15.0, 15.0, 15.0, 96, 6, 25, 12, 3, 1,   4,    3,  2,
           4,  1.5, 0.75,  0.35, 150,  1.2,  0.08, 6,  1, 3,  1,  4, 1.5, 0.45, 0.2},
          {19, 2,   1.1, 0.3, 10.0, 20.0, 30.0, 32, 1, 10, 5, 1, 0,   1,   1,  1,
           1,  1.1, 0.3, 0.3, 50,   0.5,  0.02, 2,  0, 1,  0, 1, 1.1, 0.2, 0.2},
      };

  for (const auto &params : parameters) {
    autopas::LiveInfo info(std::get<0>(params),   // numOwnedParticles
                           std::get<1>(params),   // numHaloParticles
                           std::get<2>(params),   // cutoff
                           std::get<3>(params),   // skin
                           std::get<4>(params),   // domainSizeX
                           std::get<5>(params),   // domainSizeY
                           std::get<6>(params),   // domainSizeZ
                           std::get<7>(params),   // particleSize
                           std::get<8>(params),   // threadCount
                           std::get<9>(params),   // rebuildFrequency
                           std::get<10>(params),  // numCells
                           std::get<11>(params),  // numEmptyCells
                           std::get<12>(params),  // minParticlesPerCell
                           std::get<13>(params),  // maxParticlesPerCell
                           std::get<14>(params),  // medianParticlesPerCell
                           std::get<15>(params),  // lowerQuartileParticlesPerCell
                           std::get<16>(params),  // upperQuartileParticlesPerCell
                           std::get<17>(params),  // meanParticlesPerCell
                           std::get<18>(params),  // particlesPerCellStdDev
                           std::get<19>(params),  // relativeParticlesPerCellStdDev
                           std::get<20>(params),  // estimatedNumNeighborInteractions
                           std::get<21>(params),  // particleDependentBinMaxDensity
                           std::get<22>(params),  // particleDependentBinDensityStdDev
                           std::get<23>(params),  // maxParticlesPerBlurredBin
                           std::get<24>(params),  // minParticlesPerBlurredBin
                           std::get<25>(params),  // medianParticlesPerBlurredBin
                           std::get<26>(params),  // lowerQuartileParticlesPerBlurredBin
                           std::get<27>(params),  // upperQuartileParticlesPerBlurredBin
                           std::get<28>(params),  // meanParticlesPerBlurredBin
                           std::get<29>(params),  // particlesPerBlurredBinStdDev
                           std::get<30>(params)   // relativeParticlesPerBlurredBinStdDev
    );

    EXPECT_EQ(info.get<size_t>("numOwnedParticles"), std::get<0>(params));
    EXPECT_EQ(info.get<size_t>("numHaloParticles"), std::get<1>(params));
    EXPECT_EQ(info.get<double>("cutoff"), std::get<2>(params));
    EXPECT_EQ(info.get<double>("skin"), std::get<3>(params));
    EXPECT_EQ(info.get<double>("domainSizeX"), std::get<4>(params));
    EXPECT_EQ(info.get<double>("domainSizeY"), std::get<5>(params));
    EXPECT_EQ(info.get<double>("domainSizeZ"), std::get<6>(params));
    EXPECT_EQ(info.get<size_t>("particleSize"), std::get<7>(params));
    EXPECT_EQ(info.get<size_t>("threadCount"), std::get<8>(params));
    EXPECT_EQ(info.get<size_t>("rebuildFrequency"), std::get<9>(params));
    EXPECT_EQ(info.get<size_t>("numCells"), std::get<10>(params));
    EXPECT_EQ(info.get<size_t>("numEmptyCells"), std::get<11>(params));
    EXPECT_EQ(info.get<size_t>("minParticlesPerCell"), std::get<12>(params));
    EXPECT_EQ(info.get<size_t>("maxParticlesPerCell"), std::get<13>(params));
    EXPECT_EQ(info.get<size_t>("medianParticlesPerCell"), std::get<14>(params));
    EXPECT_EQ(info.get<size_t>("lowerQuartileParticlesPerCell"), std::get<15>(params));
    EXPECT_EQ(info.get<size_t>("upperQuartileParticlesPerCell"), std::get<16>(params));
    EXPECT_EQ(info.get<double>("meanParticlesPerCell"), std::get<17>(params));
    EXPECT_EQ(info.get<double>("particlesPerCellStdDev"), std::get<18>(params));
    EXPECT_EQ(info.get<double>("relativeParticlesPerCellStdDev"), std::get<19>(params));
    EXPECT_EQ(info.get<size_t>("estimatedNumNeighborInteractions"), std::get<20>(params));
    EXPECT_EQ(info.get<double>("particleDependentBinMaxDensity"), std::get<21>(params));
    EXPECT_EQ(info.get<double>("particleDependentBinDensityStdDev"), std::get<22>(params));
    EXPECT_EQ(info.get<size_t>("maxParticlesPerBlurredBin"), std::get<23>(params));
    EXPECT_EQ(info.get<size_t>("minParticlesPerBlurredBin"), std::get<24>(params));
    EXPECT_EQ(info.get<size_t>("medianParticlesPerBlurredBin"), std::get<25>(params));
    EXPECT_EQ(info.get<size_t>("lowerQuartileParticlesPerBlurredBin"), std::get<26>(params));
    EXPECT_EQ(info.get<size_t>("upperQuartileParticlesPerBlurredBin"), std::get<27>(params));
    EXPECT_EQ(info.get<double>("meanParticlesPerBlurredBin"), std::get<28>(params));
    EXPECT_EQ(info.get<double>("particlesPerBlurredBinStdDev"), std::get<29>(params));
    EXPECT_EQ(info.get<double>("relativeParticlesPerBlurredBinStdDev"), std::get<30>(params));
  }
}

/**
 * Tests that the constructor that takes the live infos as arguments initializes all required infos.
 *
 * This test verifies that all live infos collected by the gather function are also set by the constructor which takes
 * them as arguments.
 */
TEST_F(LiveInfoTest, ConstructorNumInfos) {
  // The values here are arbitrary as we only test the number of infos collected.
  autopas::LiveInfo infoCtor(1, 0, 0.1, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 1, 0.0, 0.0, 1,
                             0, 0, 0, 0, 0.0, 0.0, 0.0);
  autopas::LiveInfo infoGather;
  autopas::AutoPas<ParticleFP64> container;
  container.setCutoff(0.5);
  container.setVerletSkin(0.4);
  container.setBoxMin({0., 0., 0.});
  container.setBoxMax({20., 20., 20.});
  container.setVerletRebuildFrequency(20);
  container.init();
  infoGather.gather(container.begin(), container.getVerletRebuildFrequency(), container.getNumberOfParticles(),
                    container.getBoxMin(), container.getBoxMax(), container.getCutoff(), container.getVerletSkin());
  EXPECT_EQ(infoCtor.get().size(), infoGather.get().size())
      << "Number of infos in constructor and gather should match.";
}

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
  const auto numCellBinsPerDim = floorAndCast<size_t>(boxSize / (cutoff + skin));
  const auto numCells = numCellBinsPerDim[0] * numCellBinsPerDim[1] * numCellBinsPerDim[2];
  const auto cellBinDims = boxSize / staticCastArray<double>(numCellBinsPerDim);
  auto cellBins = autopas::utils::ParticleBinStructure(numCellBinsPerDim, cellBinDims, boxMin, boxMax, cutoff);

  // PD = Particle Dependent
  const auto targetNumberOfPDBins = std::ceil(static_cast<double>(numParticlesOwned) / 10.);
  const auto targetNumberOfPDBinsPerDim = std::cbrt(targetNumberOfPDBins);
  const auto numberOfPDBinsPerDim = static_cast<size_t>(std::floor(targetNumberOfPDBinsPerDim));
  const auto PDBinDimensions = boxSize / static_cast<double>(numberOfPDBinsPerDim);
  auto particleDependentBins =
      autopas::utils::ParticleBinStructure(numberOfPDBinsPerDim, PDBinDimensions, boxMin, boxMax, cutoff);

  constexpr auto blurredBinDims = boxSize / std::array<double, 3>{3., 3., 3.};
  auto blurredBins = autopas::utils::ParticleBinStructure({3, 3, 3}, blurredBinDims, boxMin, boxMax, cutoff);

  // Add particles to the container and bin them
  size_t id{0};
  for (size_t i = 0; i < numParticlesOwned; ++i, ++id) {
    ParticleFP64 p({5., 5., 5.}, {0., 0., 0.}, id);
    container.addParticle(p);
    cellBins.countParticle({5., 5., 5.});
    particleDependentBins.countParticle({5., 5., 5.});
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
  liveInfo.gather(container.begin(), container.getVerletRebuildFrequency(), container.getNumberOfParticles(),
                  container.getBoxMin(), container.getBoxMax(), container.getCutoff(), container.getVerletSkin());

  cellBins.calculateStatistics();
  particleDependentBins.calculateStatistics();
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
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("estimatedNumNeighborInteractions"),
                   cellBins.getEstimatedNumberOfNeighborInteractions());

  EXPECT_EQ(liveInfo.get<double>("particleDependentBinMaxDensity"), particleDependentBins.getMaxDensity());
  EXPECT_EQ(liveInfo.get<double>("particleDependentBinDensityStdDev"), particleDependentBins.getStdDevDensity());

  EXPECT_EQ(liveInfo.get<size_t>("maxParticlesPerBlurredBin"), blurredBins.getMaxParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("minParticlesPerBlurredBin"), blurredBins.getMinParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("medianParticlesPerBlurredBin"), blurredBins.getMedianParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("lowerQuartileParticlesPerBlurredBin"), blurredBins.getLowerQuartileParticlesPerBin());
  EXPECT_EQ(liveInfo.get<size_t>("upperQuartileParticlesPerBlurredBin"), blurredBins.getUpperQuartileParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("meanParticlesPerBlurredBin"), blurredBins.getMeanParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("relativeParticlesPerBlurredBinStdDev"),
                   blurredBins.getRelStdDevParticlesPerBin());
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("particlesPerBlurredBinStdDev"), blurredBins.getStdDevParticlesPerBin());
}

/**
 * Tests that no LiveInfo Statistics are missing from the above test.
 *
 * There should be
 * - 11 directly calculated statistics
 * - 10 cellBin statistics
 * - 2 particleDependentBin statistics
 * - 8 blurredBin statistics
 *
 * = 31
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
  liveInfo.gather(container.begin(), container.getVerletRebuildFrequency(), container.getNumberOfParticles(),
                  container.getBoxMin(), container.getBoxMax(), container.getCutoff(), container.getVerletSkin());

  EXPECT_EQ(size(liveInfo.get()), 31) << "LiveInfoTest is missing a test for one of the LiveInfo statistics, or"
                                         " includes a test for a statistic that no longer exists.";
}

/**
 * Tests LiveInfo in the case that there are no particles.
 */
TEST_F(LiveInfoTest, NoParticleTest) {
  using namespace autopas::utils::ArrayMath;
  using namespace literals;

  autopas::LiveInfo liveInfo;

  constexpr double cutoff{2.0};
  constexpr double skin{0.1};
  constexpr std::array<double, 3> boxMin{0., 0., 0.};
  constexpr std::array<double, 3> boxMax{10., 10., 10.};
  constexpr auto boxSize = boxMax - boxMin;
  constexpr size_t rebuildFrequency{10};

  const auto numCellBinsPerDim = floorAndCast<size_t>(boxSize / (cutoff + skin));
  const auto numCells = numCellBinsPerDim[0] * numCellBinsPerDim[1] * numCellBinsPerDim[2];

  auto container = autopas::AutoPas<ParticleFP64>();

  container.setCutoff(cutoff);
  container.setVerletSkin(skin);
  container.setBoxMin(boxMin);
  container.setBoxMax(boxMax);
  container.setVerletRebuildFrequency(rebuildFrequency);

  container.init();

  // Gather information from the container
  liveInfo.gather(container.begin(), container.getVerletRebuildFrequency(), container.getNumberOfParticles(),
                  container.getBoxMin(), container.getBoxMax(), container.getCutoff(), container.getVerletSkin());

  // Compare directly calculated LiveInfo statistics
  EXPECT_EQ(liveInfo.get<size_t>("numOwnedParticles"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("numHaloParticles"), 0);
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
  EXPECT_EQ(liveInfo.get<size_t>("numEmptyCells"), numCells);
  EXPECT_EQ(liveInfo.get<size_t>("maxParticlesPerCell"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("minParticlesPerCell"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("medianParticlesPerCell"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("lowerQuartileParticlesPerCell"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("upperQuartileParticlesPerCell"), 0);
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("meanParticlesPerCell"), 0.);
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("relativeParticlesPerCellStdDev"), 0.);
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("particlesPerCellStdDev"), 0.);
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("estimatedNumNeighborInteractions"), 0.);

  EXPECT_EQ(liveInfo.get<double>("particleDependentBinMaxDensity"), 0);
  EXPECT_EQ(liveInfo.get<double>("particleDependentBinDensityStdDev"), 0);

  EXPECT_EQ(liveInfo.get<size_t>("maxParticlesPerBlurredBin"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("minParticlesPerBlurredBin"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("medianParticlesPerBlurredBin"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("lowerQuartileParticlesPerBlurredBin"), 0);
  EXPECT_EQ(liveInfo.get<size_t>("upperQuartileParticlesPerBlurredBin"), 0);
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("meanParticlesPerBlurredBin"), 0.);
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("relativeParticlesPerBlurredBinStdDev"), 0.);
  EXPECT_DOUBLE_EQ(liveInfo.get<double>("particlesPerBlurredBinStdDev"), 0.);
}