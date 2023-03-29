/**
 * @file AutoPasAllContainersTest.cpp
 * @author F. Gratl
 * @date 29.03.23
 */

#include "AutoPasAllContainersTest.h"

#include "autopas/AutoPasImpl.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;

INSTANTIATE_TEST_SUITE_P(Generated, AutoPasAllContainersTest,
                         ::testing::ValuesIn(autopas::ContainerOption::getAllOptions()));

/**
 * Fill a container with randomly generated particles, invoke addParticles, then check if the expected amount was added.
 */
TEST_P(AutoPasAllContainersTest, addParticlesTest) {
  /// Config
  constexpr auto numMols = 100;
  constexpr std::array<double, 3> boxMin{0., 0., 0.};
  constexpr std::array<double, 3> boxMax{10., 10., 10.};

  /// Setup
  // initialize AutoPas
  autopas::AutoPas<Molecule> autoPas;
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  const auto containerOption = GetParam();
  autoPas.setAllowedContainers({containerOption});
  autoPas.setAllowedTraversals(autopas::TraversalOption::getAllOptions());
  autoPas.setNumSamples(autoPas.getTuningInterval());  // to avoid warnings
  autoPas.init();
  // generate some random molecules
  std::vector<Molecule> mols;
  mols.reserve(numMols);
  for (size_t i = 0; i < numMols; ++i) {
    mols.emplace_back(autopasTools::generators::RandomGenerator::randomPosition(boxMin, boxMax),
                      std::array<double, 3>{0., 0., 0.}, i, 0);
  }

  /// Test
  // add every second particle
  autoPas.addParticlesIf(mols, [](const auto &p) { return p.getID() % 2 == 0; });

  /// Expectations
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned), numMols / 2)
      << "AutoPas' counters are wrong!";
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::halo), 0) << "AutoPas' counters are wrong!";
  // count what is actually in the container
  size_t numParticlesFound = 0;
  size_t numHalosFound = 0;
  for (const auto &p : autoPas) {
    if (p.isOwned()) {
      ++numParticlesFound;
    }
    if (p.isHalo()) {
      ++numHalosFound;
    }
  }
  EXPECT_EQ(numParticlesFound, numMols / 2) << "Either iterators are broken or the wrong number of particles is wrong!";
  EXPECT_EQ(numHalosFound, 0) << "Either iterators are broken or the wrong number of particles is wrong!";
}
