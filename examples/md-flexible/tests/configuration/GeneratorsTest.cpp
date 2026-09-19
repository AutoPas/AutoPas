/**
 * @file GeneratorsTest.cpp
 * @author N. Fottner
 * @date 3/8/19
 */

#include "GeneratorsTest.h"

#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/generators/GridGenerator.h"
#include "src/configuration/YamlParser.h"
#include "src/configuration/objects/CubeClosestPacked.h"
#include "src/configuration/objects/CubeGauss.h"
#include "src/configuration/objects/CubeGrid.h"
#include "src/configuration/objects/CubeUniform.h"
#include "src/configuration/objects/Sphere.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * This test checks if the GridGenerator fills the container with particles that are inside the box.
 */
TEST_F(GeneratorsTest, GridFillwithBoxMin) {
  auto autoPas = autopas::AutoPas<ParticleType>(std::cout);
  constexpr std::array<double, 3> boxMin = {5., 5., 5.};
  constexpr std::array<double, 3> boxMax = {10., 10., 10.};
  autoPas.setBoxMax(boxMax);
  autoPas.setBoxMin(boxMin);
  const ParticleType dummy;

  autoPas.init();
  autopas::generators::GridGenerator::fillWithParticles(autoPas, {5, 5, 5}, dummy, {1, 1, 1}, boxMin);
  AUTOPAS_OPENMP(parallel)
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_TRUE(autopas::utils::inBox(iter->getR(), boxMin, boxMax));
  }
}

/**
 * This test expects multipleObjectsWithMultipleTypesTest.yaml to be placed in md-flexible/tests/yamlTestFiles
 */
TEST_F(GeneratorsTest, MultipleObjectGeneration) {
  int myRank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  if (myRank != 0) {
    GTEST_SKIP() << "[Rank " << myRank
                 << "] MultipleObjectGeneration test works only on rank 0 because MDFlexConfig only generates "
                    "particles on rank 0.";
  }

  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  EXPECT_THAT(configuration.getObjectsByType<CubeGrid>(), ::testing::SizeIs(1));
  EXPECT_THAT(configuration.getObjectsByType<CubeGauss>(), ::testing::SizeIs(1));
  EXPECT_THAT(configuration.getObjectsByType<CubeUniform>(), ::testing::SizeIs(1));
  EXPECT_THAT(configuration.getObjectsByType<Sphere>(), ::testing::SizeIs(1));
  EXPECT_THAT(configuration.getObjectsByType<CubeClosestPacked>(), ::testing::SizeIs(1));

  // counters to checks if all particles types are well initialized for different Objects:
  int gridCounter = 0;
  int gaussCounter = 0;
  int uniformCounter = 0;
  int sphereCounter = 0;
  int closestCounter = 0;

  const std::array<double, 3> velocity = {0., 0., 0.};
  for (auto &particle : configuration.particles) {
    EXPECT_EQ(velocity, particle.getV());  // velocity set to {0.,0.,0.} in parsingFile
    switch (particle.getTypeId()) {
      case 0: {
        gridCounter++;
        break;
      }
      case 1: {
        gaussCounter++;
        break;
      }
      case 2: {
        uniformCounter++;
        break;
      }
      case 3: {
        sphereCounter++;
        break;
      }
      case 4: {
        closestCounter++;
        break;
      }
      default: {
        throw std::runtime_error("something went wrong with the Types");
      }
    }
  }

  EXPECT_EQ(gridCounter, configuration.getObjectsByType<CubeGrid>().at(0)->getParticlesTotal());
  EXPECT_EQ(gaussCounter, configuration.getObjectsByType<CubeGauss>().at(0)->getParticlesTotal());
  EXPECT_EQ(uniformCounter, configuration.getObjectsByType<CubeUniform>().at(0)->getParticlesTotal());
  EXPECT_EQ(sphereCounter, configuration.getObjectsByType<Sphere>().at(0)->getParticlesTotal());
  EXPECT_EQ(closestCounter, configuration.getObjectsByType<CubeClosestPacked>().at(0)->getParticlesTotal());
  // check if during initialization, not 2 Particles were initialized with same id
  std::set<size_t> ids;
  for (auto &particle : configuration.particles) {
    const auto particleId = particle.getID();
    ASSERT_EQ(ids.count(particleId), 0) << "Two particles have the same ID " << particleId;
    ids.insert(particleId);
  }
  EXPECT_EQ(ids.size(), configuration.particles.size());
}

/**
 * This test checks if the CubeClosestPacked generator with FCC lattice structure fills the container with particles
 * that are inside the box.
 */
TEST_F(GeneratorsTest, CubeClosestPackedFCC) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr double spacing = 1.0;
  constexpr std::array<double, 3> boxLength = {4.0, 4.0, 4.0};
  constexpr std::array<double, 3> bottomLeft = {0., 0., 0.};

  const CubeClosestPacked cube(velocity, typeId, spacing, boxLength, bottomLeft,
                               CubeClosestPacked::LatticeStructure::FCC);
  EXPECT_EQ(cube.getLatticeStructure(), CubeClosestPacked::LatticeStructure::FCC);

  std::vector<ParticleType> particles;
  cube.generate(particles);

  EXPECT_EQ(particles.size(), cube.getParticlesTotal());
  EXPECT_GT(particles.size(), 0);

  // Verify all particles are inside bounds
  for (const auto &p : particles) {
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_GE(p.getR()[d], bottomLeft[d]);
      EXPECT_LT(p.getR()[d], bottomLeft[d] + boxLength[d]);
    }
  }

  // Check nearest-neighbor distance among first few particles
  double minDistance = std::numeric_limits<double>::max();
  for (size_t i = 0; i < std::min<size_t>(particles.size(), 30); ++i) {
    for (size_t j = i + 1; j < std::min<size_t>(particles.size(), 30); ++j) {
      const auto diff = autopas::utils::ArrayMath::sub(particles[i].getR(), particles[j].getR());
      const double dist = autopas::utils::ArrayMath::dot(diff, diff);
      minDistance = std::min(minDistance, std::sqrt(dist));
    }
  }
  EXPECT_NEAR(minDistance, spacing, 1e-5);
}

/**
 * This test checks if the CubeClosestPacked generator with HCP lattice structure fills the container with particles
 * that are inside the box.
 */
TEST_F(GeneratorsTest, CubeClosestPackedHCP) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr double spacing = 1.0;
  constexpr std::array<double, 3> boxLength = {4.0, 4.0, 4.0};
  constexpr std::array<double, 3> bottomLeft = {0., 0., 0.};

  const CubeClosestPacked cube(velocity, typeId, spacing, boxLength, bottomLeft,
                               CubeClosestPacked::LatticeStructure::HCP);
  EXPECT_EQ(cube.getLatticeStructure(), CubeClosestPacked::LatticeStructure::HCP);

  std::vector<ParticleType> particles;
  cube.generate(particles);

  EXPECT_EQ(particles.size(), cube.getParticlesTotal());
  EXPECT_GT(particles.size(), 0);

  for (const auto &p : particles) {
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_GE(p.getR()[d], bottomLeft[d]);
      EXPECT_LT(p.getR()[d], bottomLeft[d] + boxLength[d]);
    }
  }
}

/**
 * Tests CubeGrid generator across multiple scenarios of densities, dimensions, alignments, and origins.
 */
TEST_F(GeneratorsTest, CubeGridDensityScenarios) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;

  const std::vector<GridDensityScenario> gridScenarios = {
      {{2, 2, 2}, {0.0, 0.0, 0.0}, 8.0, true, "Grid cubic 2x2x2 centered"},
      {{2, 2, 2}, {0.0, 0.0, 0.0}, 8.0, false, "Grid cubic 2x2x2 uncentered"},
      {{3, 3, 3}, {1.0, 2.0, 3.0}, 1.0, true, "Grid cubic unit density centered"},
      {{3, 3, 3}, {1.0, 2.0, 3.0}, 1.0, false, "Grid cubic unit density uncentered"},
      {{2, 2, 2}, {1.0, 2.0, 3.0}, 0.125, true, "Grid cubic low density centered"},
      {{4, 2, 3}, {1.0, -2.0, 3.0}, 1.0, true, "Grid anisotropic centered"},
      {{4, 2, 3}, {1.0, -2.0, 3.0}, 1.0, false, "Grid anisotropic uncentered"},
      {{5, 3, 2}, {-2.0, 1.5, -0.5}, 2.0, true, "Grid anisotropic negative origin centered"},
      {{5, 3, 2}, {-2.0, 1.5, -0.5}, 2.0, false, "Grid anisotropic negative origin uncentered"},
  };

  for (const auto &[particlesPerDim, bottomLeft, density, centered, description] : gridScenarios) {
    SCOPED_TRACE(description);
    const CubeGrid grid(velocity, typeId, particlesPerDim, bottomLeft, density, centered);
    const double expectedSpacing = std::cbrt(1.0 / density);
    const size_t expectedTotal = particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2];

    EXPECT_DOUBLE_EQ(grid.getParticleDensity(), density);
    EXPECT_DOUBLE_EQ(grid.getParticleSpacing(), expectedSpacing);
    EXPECT_EQ(grid.getParticlesTotal(), expectedTotal);

    const auto boxMin = grid.getBoxMin();
    const auto boxMax = grid.getBoxMax();
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_DOUBLE_EQ(boxMin[d], bottomLeft[d]);
      const double span =
          static_cast<double>(centered ? particlesPerDim[d] : (particlesPerDim[d] - 1)) * expectedSpacing;
      EXPECT_NEAR(boxMax[d], bottomLeft[d] + span, 1e-10);
    }

    if (centered) {
      const double boxVolume = (boxMax[0] - boxMin[0]) * (boxMax[1] - boxMin[1]) * (boxMax[2] - boxMin[2]);
      EXPECT_NEAR(static_cast<double>(grid.getParticlesTotal()) / boxVolume, density, 1e-10);
    }

    std::vector<ParticleType> particles;
    grid.generate(particles);
    ASSERT_EQ(particles.size(), expectedTotal);

    // Verify alignment of first particle
    for (size_t d = 0; d < 3; ++d) {
      const double expectedCoord = bottomLeft[d] + (centered ? 0.5 * expectedSpacing : 0.0);
      EXPECT_NEAR(particles[0].getR()[d], expectedCoord, 1e-10);
    }

    // Verify bounds
    for (const auto &p : particles) {
      for (size_t d = 0; d < 3; ++d) {
        EXPECT_GE(p.getR()[d], boxMin[d]);
        EXPECT_LE(p.getR()[d], boxMax[d]);
      }
    }
  }
}

/**
 * Tests CubeUniform generator across multiple scenarios of densities, dimensions, roundings, and origins.
 */
TEST_F(GeneratorsTest, CubeUniformDensityScenarios) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;

  const std::vector<UniformDensityScenario> uniformScenarios = {
      {{2.0, 2.0, 2.0}, {0.0, 0.0, 0.0}, 2.5, "Uniform cubic exact integer"},
      {{4.0, 4.0, 4.0}, {0.0, 0.0, 0.0}, 1.0, "Uniform cubic unit density"},
      {{4.0, 4.0, 4.0}, {1.0, 2.0, 3.0}, 1.31, "Uniform cubic round down"},
      {{4.0, 4.0, 4.0}, {1.0, 2.0, 3.0}, 1.36, "Uniform cubic round up"},
      {{3.5, 5.0, 4.0}, {-1.0, 2.0, -0.5}, 0.8, "Uniform anisotropic low density"},
      {{3.5, 5.0, 4.0}, {-1.0, 2.0, -0.5}, 1.5, "Uniform anisotropic high density"},
  };

  for (const auto &[boxLength, bottomLeft, density, description] : uniformScenarios) {
    SCOPED_TRACE(description);
    const CubeUniform uniform(velocity, typeId, boxLength, bottomLeft, density);
    const double volume = boxLength[0] * boxLength[1] * boxLength[2];
    const size_t expectedParticles = static_cast<size_t>(std::round(density * volume));
    const double expectedActualDensity = static_cast<double>(expectedParticles) / volume;

    EXPECT_EQ(uniform.getParticlesTotal(), expectedParticles);
    EXPECT_DOUBLE_EQ(uniform.getParticleDensity(), expectedActualDensity);

    std::vector<ParticleType> particles;
    uniform.generate(particles);
    ASSERT_EQ(particles.size(), expectedParticles);

    for (const auto &p : particles) {
      for (size_t d = 0; d < 3; ++d) {
        EXPECT_GE(p.getR()[d], bottomLeft[d]);
        EXPECT_LT(p.getR()[d], bottomLeft[d] + boxLength[d]);
      }
    }
  }
}

/**
 * Tests CubeClosestPacked generator across multiple scenarios of densities, dimensions, alignments, and lattice
 * structures.
 */
TEST_F(GeneratorsTest, CubeClosestPackedDensityScenarios) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;

  const std::vector<ClosestPackedDensityScenario> scenarios = {
      // FCC Cubic boxes across diverse densities and alignments
      {CubeClosestPacked::LatticeStructure::FCC,
       {4.0, 4.0, 4.0},
       {0.0, 0.0, 0.0},
       0.5,
       true,
       "FCC cubic low density centered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {4.0, 4.0, 4.0},
       {1.0, 2.0, 3.0},
       0.984375,
       true,
       "FCC cubic exact density centered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {4.0, 4.0, 4.0},
       {1.0, 2.0, 3.0},
       0.984375,
       false,
       "FCC cubic exact density uncentered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {4.0, 4.0, 4.0},
       {0.0, 0.0, 0.0},
       1.0,
       true,
       "FCC cubic unit density centered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {4.0, 4.0, 4.0},
       {0.0, 0.0, 0.0},
       1.0,
       false,
       "FCC cubic unit density uncentered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {4.0, 4.0, 4.0},
       {-1.0, 1.0, -2.0},
       2.0,
       true,
       "FCC cubic high density centered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {4.0, 4.0, 4.0},
       {-1.0, 1.0, -2.0},
       2.0,
       false,
       "FCC cubic high density uncentered"},

      // FCC Anisotropic boxes
      {CubeClosestPacked::LatticeStructure::FCC,
       {3.5, 5.0, 4.2},
       {1.0, -2.0, 3.0},
       1.5,
       true,
       "FCC anisotropic centered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {3.5, 5.0, 4.2},
       {1.0, -2.0, 3.0},
       1.5,
       false,
       "FCC anisotropic uncentered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {5.0, 2.0, 3.0},
       {0.0, 0.0, 0.0},
       0.8,
       true,
       "FCC anisotropic low density centered"},
      {CubeClosestPacked::LatticeStructure::FCC,
       {2.5, 4.0, 3.2},
       {-2.0, 0.5, 1.5},
       2.5,
       false,
       "FCC anisotropic high density uncentered"},

      // HCP Cubic boxes across diverse densities and alignments
      {CubeClosestPacked::LatticeStructure::HCP,
       {4.0, 4.0, 4.0},
       {0.0, 0.0, 0.0},
       0.5,
       true,
       "HCP cubic low density centered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {4.0, 4.0, 4.0},
       {1.0, 2.0, 3.0},
       1.0,
       true,
       "HCP cubic unit density centered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {4.0, 4.0, 4.0},
       {1.0, 2.0, 3.0},
       1.0,
       false,
       "HCP cubic unit density uncentered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {4.0, 4.0, 4.0},
       {0.0, 0.0, 0.0},
       1.09375,
       true,
       "HCP cubic exact density centered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {4.0, 4.0, 4.0},
       {-1.0, 1.0, -2.0},
       2.0,
       true,
       "HCP cubic high density centered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {4.0, 4.0, 4.0},
       {-1.0, 1.0, -2.0},
       2.0,
       false,
       "HCP cubic high density uncentered"},

      // HCP Anisotropic boxes
      {CubeClosestPacked::LatticeStructure::HCP,
       {3.5, 5.0, 4.2},
       {1.0, -2.0, 3.0},
       1.5,
       true,
       "HCP anisotropic centered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {3.5, 5.0, 4.2},
       {1.0, -2.0, 3.0},
       1.5,
       false,
       "HCP anisotropic uncentered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {5.0, 2.0, 3.0},
       {0.0, 0.0, 0.0},
       0.8,
       true,
       "HCP anisotropic low density centered"},
      {CubeClosestPacked::LatticeStructure::HCP,
       {2.5, 4.0, 3.2},
       {-2.0, 0.5, 1.5},
       2.5,
       false,
       "HCP anisotropic high density uncentered"},
  };

  for (const auto &[structure, boxLength, bottomLeft, density, centered, description] : scenarios) {
    SCOPED_TRACE(description);
    const CubeClosestPacked cube(velocity, typeId, boxLength, bottomLeft, density, structure, centered);

    const double s0 = std::cbrt(std::sqrt(2.0) / density);
    const double volume = boxLength[0] * boxLength[1] * boxLength[2];

    EXPECT_EQ(cube.getLatticeStructure(), structure);
    // Guarantee no tighter spacing: s* >= s0
    EXPECT_GE(cube.getParticleSpacing(), s0 - 1e-12);
    EXPECT_GT(cube.getParticlesTotal(), 0);
    EXPECT_DOUBLE_EQ(cube.getParticleDensity(), static_cast<double>(cube.getParticlesTotal()) / volume);

    std::vector<ParticleType> particles;
    cube.generate(particles);
    ASSERT_EQ(particles.size(), cube.getParticlesTotal());

    // Verify all particles inside bounding box [bottomLeft, bottomLeft + boxLength)
    for (const auto &p : particles) {
      for (size_t d = 0; d < 3; ++d) {
        EXPECT_GE(p.getR()[d], bottomLeft[d]);
        EXPECT_LT(p.getR()[d], bottomLeft[d] + boxLength[d]);
      }
    }

    // Verify alignment offset of first particle
    const double spacing = cube.getParticleSpacing();
    if (not centered) {
      for (size_t d = 0; d < 3; ++d) {
        EXPECT_NEAR(particles[0].getR()[d], bottomLeft[d], 1e-10);
      }
    } else {
      if (structure == CubeClosestPacked::LatticeStructure::FCC) {
        const double a = std::sqrt(2.0) * spacing;
        for (size_t d = 0; d < 3; ++d) {
          EXPECT_NEAR(particles[0].getR()[d], bottomLeft[d] + a / 4.0, 1e-10);
        }
      } else {
        const double yOffset = spacing * std::sqrt(1. / 12.);
        const double spacingLayer = spacing * std::sqrt(2. / 3.);
        EXPECT_NEAR(particles[0].getR()[0], bottomLeft[0] + spacing / 4.0, 1e-10);
        EXPECT_NEAR(particles[0].getR()[1], bottomLeft[1] + yOffset, 1e-10);
        EXPECT_NEAR(particles[0].getR()[2], bottomLeft[2] + spacingLayer / 2.0, 1e-10);
      }
    }

    // Verify nearest-neighbor distance >= spacing - 1e-5
    if (particles.size() > 1) {
      double minDistance = std::numeric_limits<double>::max();
      const size_t sampleCount = std::min<size_t>(particles.size(), 30);
      for (size_t i = 0; i < sampleCount; ++i) {
        for (size_t j = i + 1; j < sampleCount; ++j) {
          const auto diff = autopas::utils::ArrayMath::sub(particles[i].getR(), particles[j].getR());
          const double dist = autopas::utils::ArrayMath::dot(diff, diff);
          minDistance = std::min(minDistance, std::sqrt(dist));
        }
      }
      EXPECT_NEAR(minDistance, spacing, 1e-5);
      EXPECT_GE(minDistance, s0 - 1e-5);
    }
  }
}

/**
 * This test checks if particle IDs are continuous and unique across multiple generators constructed by density.
 */
TEST_F(GeneratorsTest, DensityIDContinuityAndCumulativeGeneration) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;

  std::vector<ParticleType> particles;

  // 1. Generate CubeGrid with density
  const CubeGrid grid(velocity, typeId, {2, 2, 2}, {0.0, 0.0, 0.0}, 2.0);
  const size_t gridCount = grid.getParticlesTotal();
  grid.generate(particles);
  EXPECT_EQ(particles.size(), gridCount);

  // 2. Generate CubeUniform with density
  const CubeUniform uniform(velocity, typeId, {2.0, 2.0, 2.0}, {10.0, 0.0, 0.0}, 1.5);
  const size_t uniformCount = uniform.getParticlesTotal();
  uniform.generate(particles);
  EXPECT_EQ(particles.size(), gridCount + uniformCount);

  // 3. Generate CubeClosestPacked (FCC) with density
  const CubeClosestPacked ccpFCC(velocity, typeId, {3.0, 3.0, 3.0}, {20.0, 0.0, 0.0}, 1.0,
                                 CubeClosestPacked::LatticeStructure::FCC);
  const size_t fccCount = ccpFCC.getParticlesTotal();
  ccpFCC.generate(particles);
  EXPECT_EQ(particles.size(), gridCount + uniformCount + fccCount);

  // 4. Generate CubeClosestPacked (HCP) with density
  const CubeClosestPacked ccpHCP(velocity, typeId, {3.0, 3.0, 3.0}, {30.0, 0.0, 0.0}, 1.0,
                                 CubeClosestPacked::LatticeStructure::HCP);
  const size_t hcpCount = ccpHCP.getParticlesTotal();
  ccpHCP.generate(particles);
  EXPECT_EQ(particles.size(), gridCount + uniformCount + fccCount + hcpCount);

  // Check that IDs are strictly continuous: 0, 1, 2, ..., N - 1
  for (size_t i = 0; i < particles.size(); ++i) {
    EXPECT_EQ(particles[i].getID(), i);
  }
}

/**
 * This test checks if the spacing optimization for FCC closest packing correctly reduces density error
 * without compressing particles (s* >= s0).
 */
TEST_F(GeneratorsTest, CubeClosestPackedDensityOptimizationFCC) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr double targetDensity = 0.984375;
  constexpr std::array<double, 3> boxLength = {4.0, 4.0, 4.0};
  constexpr std::array<double, 3> bottomLeft = {0., 0., 0.};
  const double s0 = std::cbrt(std::sqrt(2.0) / targetDensity);
  constexpr size_t targetParticles = 63;  // 0.984375 * 64

  // 1. Spacing constructor with s0: unoptimized crystal packing gives 108 particles due to boundary clipping
  const CubeClosestPacked fccBySpacing(velocity, typeId, s0, boxLength, bottomLeft,
                                       CubeClosestPacked::LatticeStructure::FCC, false);
  EXPECT_DOUBLE_EQ(fccBySpacing.getParticleSpacing(), s0);
  EXPECT_EQ(fccBySpacing.getParticlesTotal(), 108);

  // 2. Density constructor: automatically optimizes spacing (s* >= s0) to achieve target density (63 particles)
  const CubeClosestPacked fccByDensity(velocity, typeId, boxLength, bottomLeft, targetDensity,
                                       CubeClosestPacked::LatticeStructure::FCC, false);
  // Guarantee NO compression: s* >= s0
  EXPECT_GE(fccByDensity.getParticleSpacing(), s0 - 1e-12);

  // Density matches target (in fact exactly 63 particles!)
  EXPECT_EQ(fccByDensity.getParticlesTotal(), targetParticles);
  EXPECT_DOUBLE_EQ(fccByDensity.getParticleDensity(), targetDensity);

  std::vector<ParticleType> particles;
  fccByDensity.generate(particles);
  EXPECT_EQ(particles.size(), targetParticles);

  // Check all particles inside box
  for (const auto &p : particles) {
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_GE(p.getR()[d], bottomLeft[d]);
      EXPECT_LT(p.getR()[d], bottomLeft[d] + boxLength[d]);
    }
  }

  // Check nearest-neighbor distance matches optimized spacing s* >= s0
  double minDistance = std::numeric_limits<double>::max();
  for (size_t i = 0; i < std::min<size_t>(particles.size(), 30); ++i) {
    for (size_t j = i + 1; j < std::min<size_t>(particles.size(), 30); ++j) {
      const auto diff = autopas::utils::ArrayMath::sub(particles[i].getR(), particles[j].getR());
      const double dist = autopas::utils::ArrayMath::dot(diff, diff);
      minDistance = std::min(minDistance, std::sqrt(dist));
    }
  }
  EXPECT_NEAR(minDistance, fccByDensity.getParticleSpacing(), 1e-5);
  EXPECT_GE(minDistance, s0 - 1e-10);
}

/**
 * This test checks if the spacing optimization for HCP closest packing correctly reduces density error
 * without compressing particles (s* >= s0).
 */
TEST_F(GeneratorsTest, CubeClosestPackedDensityOptimizationHCP) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr double targetDensity = 1.0;
  constexpr std::array<double, 3> boxLength = {4.0, 4.0, 4.0};
  constexpr std::array<double, 3> bottomLeft = {1.0, 2.0, 3.0};
  const double s0 = std::cbrt(std::sqrt(2.0) / targetDensity);
  constexpr size_t targetParticles = 64;

  // 1. Spacing constructor with s0: unoptimized packing gives 92 particles (overshoot of 28 particles)
  const CubeClosestPacked hcpBySpacing(velocity, typeId, s0, boxLength, bottomLeft,
                                       CubeClosestPacked::LatticeStructure::HCP, false);
  EXPECT_DOUBLE_EQ(hcpBySpacing.getParticleSpacing(), s0);
  const size_t unoptimizedCount = hcpBySpacing.getParticlesTotal();
  const double unoptimizedDiff = std::abs(static_cast<double>(unoptimizedCount) - targetParticles);

  // 2. Density constructor: spacing is automatically optimized (s* >= s0) to reduce error
  const CubeClosestPacked hcpByDensity(velocity, typeId, boxLength, bottomLeft, targetDensity,
                                       CubeClosestPacked::LatticeStructure::HCP, false);
  EXPECT_GE(hcpByDensity.getParticleSpacing(), s0 - 1e-12);
  const size_t optimizedCount = hcpByDensity.getParticlesTotal();
  const double optimizedDiff = std::abs(static_cast<double>(optimizedCount) - targetParticles);

  EXPECT_LE(optimizedDiff, unoptimizedDiff);

  std::vector<ParticleType> particles;
  hcpByDensity.generate(particles);
  EXPECT_EQ(particles.size(), optimizedCount);

  for (const auto &p : particles) {
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_GE(p.getR()[d], bottomLeft[d]);
      EXPECT_LT(p.getR()[d], bottomLeft[d] + boxLength[d]);
    }
  }
}

/**
 * This test checks that spacing is never decreased (particles are never compressed) when particle count
 * at s0 is already less than or equal to target count.
 */
TEST_F(GeneratorsTest, CubeClosestPackedDensityOptimizationNoCompression) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr double targetDensity = 1.0;
  constexpr std::array<double, 3> boxLength = {4.0, 4.0, 4.0};
  constexpr std::array<double, 3> bottomLeft = {0., 0., 0.};
  const double s0 = std::cbrt(std::sqrt(2.0) / targetDensity);

  // For centered HCP in 4x4x4, count at s0 is 56 <= targetParticles 64.
  // Optimization must NOT compress particles (must not reduce s below s0).
  const CubeClosestPacked hcp(velocity, typeId, boxLength, bottomLeft, targetDensity,
                              CubeClosestPacked::LatticeStructure::HCP, true);
  EXPECT_DOUBLE_EQ(hcp.getParticleSpacing(), s0);
  EXPECT_EQ(hcp.getParticlesTotal(), 56);
}

/**
 * This test checks if the CubeClosestPacked generator correctly generates particles with centered and origin alignment.
 */
TEST_F(GeneratorsTest, CubeClosestPackedAlignment) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr double spacing = 1.0;
  constexpr std::array<double, 3> bottomLeft = {1.0, 2.0, 3.0};
  constexpr std::array<double, 3> boxLength = {4.0, 4.0, 4.0};

  // 1. FCC uncentered: first particle is at bottomLeft
  {
    const CubeClosestPacked fccUncentered(velocity, typeId, spacing, boxLength, bottomLeft,
                                          CubeClosestPacked::LatticeStructure::FCC, false);
    std::vector<ParticleType> particles;
    fccUncentered.generate(particles);
    ASSERT_GT(particles.size(), 0);
    EXPECT_EQ(particles.size(), fccUncentered.getParticlesTotal());
    EXPECT_NEAR(particles[0].getR()[0], bottomLeft[0], 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], bottomLeft[1], 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], bottomLeft[2], 1e-10);
  }

  // 2. FCC centered: standoff is a / 4.0 in each dimension
  {
    const CubeClosestPacked fccCentered(velocity, typeId, spacing, boxLength, bottomLeft,
                                        CubeClosestPacked::LatticeStructure::FCC, true);
    std::vector<ParticleType> particles;
    fccCentered.generate(particles);
    ASSERT_GT(particles.size(), 0);
    EXPECT_EQ(particles.size(), fccCentered.getParticlesTotal());
    const double a = std::sqrt(2.0) * spacing;
    EXPECT_NEAR(particles[0].getR()[0], bottomLeft[0] + a / 4.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], bottomLeft[1] + a / 4.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], bottomLeft[2] + a / 4.0, 1e-10);
  }

  // 3. HCP uncentered: first particle is at bottomLeft
  {
    const CubeClosestPacked hcpUncentered(velocity, typeId, spacing, boxLength, bottomLeft,
                                          CubeClosestPacked::LatticeStructure::HCP, false);
    std::vector<ParticleType> particles;
    hcpUncentered.generate(particles);
    ASSERT_GT(particles.size(), 0);
    EXPECT_EQ(particles.size(), hcpUncentered.getParticlesTotal());
    EXPECT_NEAR(particles[0].getR()[0], bottomLeft[0], 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], bottomLeft[1], 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], bottomLeft[2], 1e-10);
  }

  // 4. HCP centered: standoff is (s / 4.0, yOffset, spacingLayer / 2.0)
  {
    const CubeClosestPacked hcpCentered(velocity, typeId, spacing, boxLength, bottomLeft,
                                        CubeClosestPacked::LatticeStructure::HCP, true);
    std::vector<ParticleType> particles;
    hcpCentered.generate(particles);
    ASSERT_GT(particles.size(), 0);
    EXPECT_EQ(particles.size(), hcpCentered.getParticlesTotal());
    const double yOffset = spacing * std::sqrt(1. / 12.);
    const double spacingLayer = spacing * std::sqrt(2. / 3.);
    EXPECT_NEAR(particles[0].getR()[0], bottomLeft[0] + spacing / 4.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], bottomLeft[1] + yOffset, 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], bottomLeft[2] + spacingLayer / 2.0, 1e-10);
  }
}

/**
 * This test checks if the CubeGrid generator correctly generates particles with centered and origin alignment.
 */
TEST_F(GeneratorsTest, CubeGridAlignmentAndGenerate) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr std::array<size_t, 3> particlesPerDim = {3, 3, 3};
  constexpr double spacing = 1.2;
  constexpr std::array<double, 3> bottomLeft = {1.0, 2.0, 3.0};

  // 1. Uncentered: first particle at bottomLeft
  {
    const CubeGrid gridUncentered(velocity, typeId, particlesPerDim, spacing, bottomLeft, false);
    std::vector<ParticleType> particles;
    gridUncentered.generate(particles);
    ASSERT_EQ(particles.size(), 27);
    EXPECT_EQ(particles.size(), gridUncentered.getParticlesTotal());
    EXPECT_NEAR(particles[0].getR()[0], bottomLeft[0], 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], bottomLeft[1], 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], bottomLeft[2], 1e-10);
  }

  // 2. Centered: first particle at bottomLeft + 0.5 * spacing
  {
    const CubeGrid gridCentered(velocity, typeId, particlesPerDim, spacing, bottomLeft, true);
    std::vector<ParticleType> particles;
    gridCentered.generate(particles);
    ASSERT_EQ(particles.size(), 27);
    EXPECT_EQ(particles.size(), gridCentered.getParticlesTotal());
    EXPECT_NEAR(particles[0].getR()[0], bottomLeft[0] + 0.5 * spacing, 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], bottomLeft[1] + 0.5 * spacing, 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], bottomLeft[2] + 0.5 * spacing, 1e-10);

    // Verify periodic distance across a periodic box: box length = particlesPerDim * spacing
    constexpr std::array<double, 3> boxLength = {3 * spacing, 3 * spacing, 3 * spacing};
    for (size_t i = 0; i < particles.size(); ++i) {
      for (size_t j = i + 1; j < particles.size(); ++j) {
        std::array<double, 3> diff = autopas::utils::ArrayMath::sub(particles[i].getR(), particles[j].getR());
        for (size_t d = 0; d < 3; ++d) {
          diff[d] -= boxLength[d] * std::round(diff[d] / boxLength[d]);
        }
        const double dist = std::sqrt(autopas::utils::ArrayMath::dot(diff, diff));
        EXPECT_GE(dist, spacing - 1e-10);
      }
    }
  }
}

/**
 * This test checks if the Sphere generator correctly generates particles within the sphere's bounding box and radial
 * cutoff.
 */
TEST_F(GeneratorsTest, SphereGenerate) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr std::array<double, 3> center = {5.0, 5.0, 5.0};
  constexpr int radius = 3;
  constexpr double spacing = 1.0;

  const Sphere sphere(velocity, typeId, center, radius, spacing);
  EXPECT_EQ(sphere.getRadius(), radius);
  EXPECT_DOUBLE_EQ(sphere.getParticleSpacing(), spacing);
  EXPECT_EQ(sphere.getCenter(), center);

  const auto boxMin = sphere.getBoxMin();
  const auto boxMax = sphere.getBoxMax();
  EXPECT_NEAR(boxMin[0], center[0] - radius * spacing, 1e-10);
  EXPECT_NEAR(boxMax[0], center[0] + radius * spacing, 1e-10);

  std::vector<ParticleType> particles;
  sphere.generate(particles);
  EXPECT_EQ(particles.size(), sphere.getParticlesTotal());
  EXPECT_GT(particles.size(), 0);

  // Check that every particle is within the sphere's bounding box and radial cutoff
  constexpr double maxRadius = (radius + 1) * spacing;
  for (const auto &p : particles) {
    const auto diff = autopas::utils::ArrayMath::sub(p.getR(), center);
    const double dist = std::sqrt(autopas::utils::ArrayMath::dot(diff, diff));
    EXPECT_LE(dist, maxRadius + 1e-10);

    for (size_t d = 0; d < 3; ++d) {
      EXPECT_GE(p.getR()[d], boxMin[d] - 1e-10);
      EXPECT_LE(p.getR()[d], boxMax[d] + 1e-10);
    }
  }
}

/**
 * This test checks if the CubeGauss generator correctly generates particles with a Gaussian distribution within the
 * specified box.
 */
TEST_F(GeneratorsTest, CubeGaussGenerate) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;
  constexpr size_t numParticles = 50;
  constexpr std::array<double, 3> boxLength = {10.0, 10.0, 10.0};
  constexpr std::array<double, 3> mean = {5.0, 5.0, 5.0};
  constexpr std::array<double, 3> stdDev = {1.0, 1.0, 1.0};
  constexpr std::array<double, 3> bottomLeft = {0.0, 0.0, 0.0};

  const CubeGauss gauss(velocity, typeId, numParticles, boxLength, mean, stdDev, bottomLeft);
  EXPECT_EQ(gauss.getParticlesTotal(), numParticles);
  EXPECT_EQ(gauss.getDistributionMean(), mean);
  EXPECT_EQ(gauss.getDistributionStdDev(), stdDev);

  std::vector<ParticleType> particles;
  gauss.generate(particles);
  EXPECT_EQ(particles.size(), numParticles);

  for (const auto &p : particles) {
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_GE(p.getR()[d], bottomLeft[d]);
      EXPECT_LT(p.getR()[d], bottomLeft[d] + boxLength[d]);
    }
  }
}

/**
 * This test checks if the particle IDs are continuous and unique across multiple generators.
 */
TEST_F(GeneratorsTest, IDContinuity) {
  constexpr std::array<double, 3> velocity = {0., 0., 0.};
  constexpr unsigned long typeId = 0;

  std::vector<ParticleType> particles;

  // Generate a CubeGrid
  const CubeGrid grid(velocity, typeId, {2, 2, 2}, 1.0, {0.0, 0.0, 0.0});
  grid.generate(particles);
  EXPECT_EQ(particles.size(), 8);

  // Generate a CubeClosestPacked into the same vector
  const CubeClosestPacked ccp(velocity, typeId, 1.0, {3.0, 3.0, 3.0}, {10.0, 10.0, 10.0});
  const size_t ccpCount = ccp.getParticlesTotal();
  ccp.generate(particles);
  EXPECT_EQ(particles.size(), 8 + ccpCount);

  // Generate a Sphere into the same vector
  const Sphere sphere(velocity, typeId, {30.0, 30.0, 30.0}, 2, 1.0);
  const size_t sphereCount = sphere.getParticlesTotal();
  sphere.generate(particles);
  EXPECT_EQ(particles.size(), 8 + ccpCount + sphereCount);

  // Check that IDs are strictly 0, 1, 2, ..., N - 1
  for (size_t i = 0; i < particles.size(); ++i) {
    EXPECT_EQ(particles[i].getID(), i);
  }
}

/**
 * This test checks if the HCP generator's particle count matches the previous implementation of getParticlesTotal()
 * when centeredAlignment was false. This ensures that the new implementation is consistent with the legacy behavior.
 */
TEST_F(GeneratorsTest, testRegressionPreviousHCPImplementation) {
  // Previous implementation of getParticlesTotal() when centeredAlignment was false
  auto legacyGetParticlesTotal = [](const std::array<double, 3> &boxLength, const double particleSpacing) -> size_t {
    const double xOffset = particleSpacing * 0.5;
    const size_t xNumRow = std::ceil(boxLength[0] / particleSpacing);
    const bool xOdd = static_cast<int>(std::ceil(boxLength[0] / xOffset)) % 2 == 1;

    const auto spacingLayer = particleSpacing * std::sqrt(2. / 3.);
    const auto spacingRow = particleSpacing * std::sqrt(3. / 4.);

    const size_t yNumEven = std::ceil(boxLength[1] / spacingRow);
    const auto yOffset = particleSpacing * std::sqrt(1. / 12.);
    const size_t yNumOdd = std::ceil((boxLength[1] - yOffset) / spacingRow);

    const size_t evenLayer = xNumRow * yNumEven - std::floor(xOdd * yNumEven * 0.5);
    const size_t oddLayer = xNumRow * yNumOdd - std::ceil(xOdd * yNumOdd * 0.5);

    const double numLayers = std::ceil(boxLength[2] / spacingLayer);
    return evenLayer * std::ceil(numLayers / 2.) + oddLayer * std::floor(numLayers / 2.);
  };

  for (const double spacing : {0.5, 1.0, 1.25, 2.3}) {
    const double spacingRow = spacing * std::sqrt(3. / 4.);
    const double spacingLayer = spacing * std::sqrt(2. / 3.);

    // Test a variety of integer and non-integer grid multiples
    for (const double fx : {0.0, 0.5, 1.0, 1.3, 2.0, 3.7, 4.0}) {
      for (const double fy : {0.0, 0.5, 1.0, 1.5, 2.0, 3.2, 5.0}) {
        for (const double fz : {0.0, 0.5, 1.0, 1.8, 2.0, 3.4, 4.0}) {
          const std::array<double, 3> boxLength = {fx * spacing, fy * spacingRow, fz * spacingLayer};
          const size_t expectedLegacy = (boxLength[0] <= 0.0 or boxLength[1] <= 0.0 or boxLength[2] <= 0.0)
                                            ? 0
                                            : legacyGetParticlesTotal(boxLength, spacing);

          const std::array<double, 3> boxMin = {0.0, 0.0, 0.0};
          const size_t actual =
              autopas::generators::HCPGenerator::getNumberOfParticles(boxMin, boxLength, spacing, false);
          EXPECT_EQ(actual, expectedLegacy) << "Mismatch with spacing=" << spacing << ", boxLength=[" << boxLength[0]
                                            << ", " << boxLength[1] << ", " << boxLength[2] << "]";
        }
      }
    }
  }
}