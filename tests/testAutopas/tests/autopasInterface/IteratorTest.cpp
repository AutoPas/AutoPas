/**
 * @file IteratorTest.cpp
 * @author seckler
 * @date 22.07.19
 */

#include "IteratorTest.h"
#include "testingHelpers/commonTypedefs.h"

constexpr double cutoff = 1.;
constexpr double skin = 0.2;
constexpr std::array<double, 3> boxMin{0., 0., 0.};
constexpr std::array<double, 3> haloBoxMin{0. - skin - cutoff, 0. - skin - cutoff, 0. - skin - cutoff};
constexpr std::array<double, 3> boxMax{10., 10., 10.};
constexpr std::array<double, 3> haloBoxMax{10. + skin + cutoff, 10. + skin + cutoff, 10. + skin + cutoff};

template <typename AutoPasT>
void defaultInit(AutoPasT &autoPas) {
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkin(skin);
  autoPas.setVerletRebuildFrequency(2);
  autoPas.setNumSamples(2);
  // init autopas
  autoPas.init();
}

/**
 * Tests the addition and iteration over particles.
 * @param containerOption
 */
void testAdditionAndIteration(testingTuple options) {
  // create AutoPas object
  autopas::AutoPas<Molecule, FMCell> autoPas;

  auto containerOption = std::get<0>(options);

  auto cellSizeOption = std::get<1>(options);

  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeOption})));

  defaultInit(autoPas);

  auto getPossible1DPositions = [&](double min, double max) -> auto {
    // ensure that all particles are at most skin away from halo!
    return std::array<double, 6>{min - cutoff, min, min + skin, max - skin, max, max + cutoff - 1e-3};
  };
  size_t id = 0;
  for (auto x : getPossible1DPositions(boxMin[0], boxMax[0])) {
    for (auto y : getPossible1DPositions(boxMin[1], boxMax[1])) {
      for (auto z : getPossible1DPositions(boxMin[2], boxMax[2])) {
        std::array<double, 3> pos{x, y, z};
        Molecule p(pos, {0., 0., 0.}, id);
        ++id;
        // add the two particles!
        if (autopas::utils::inBox(pos, boxMin, boxMax)) {
          autoPas.addParticle(p);
        } else {
          autoPas.addOrUpdateHaloParticle(p);
        }
      }
    }
  }
  {
    size_t count = 0;
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      ++count;
      EXPECT_TRUE(iter->isOwned());
    }

    EXPECT_EQ(count, 3 * 3 * 3);
  }

  // check number of halo particles
  {
    size_t count = 0;
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
      ++count;
      EXPECT_FALSE(iter->isOwned());
    }

    EXPECT_EQ(count, 6 * 6 * 6 - 3 * 3 * 3);
  }

  // check number of particles
  {
    size_t count = 0;
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
      ++count;
    }
    EXPECT_EQ(count, 6 * 6 * 6);
  }

  // check number of halo particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPas.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::haloOnly);
         iter.isValid(); ++iter) {
      ++count;
      EXPECT_FALSE(iter->isOwned());
    }

    EXPECT_EQ(count, 6 * 6 * 6 - 3 * 3 * 3);
  }

  // check number of particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPas.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::haloAndOwned);
         iter.isValid(); ++iter) {
      ++count;
    }
    EXPECT_EQ(count, 6 * 6 * 6);
  }

  // check number of owned particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPas.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::ownedOnly);
         iter.isValid(); ++iter) {
      ++count;
    }

    EXPECT_EQ(count, 3 * 3 * 3);
  }
}

TEST_P(IteratorTest, ParticleAdditionAndIteratorTestNormal) {
  auto options = GetParam();
  testAdditionAndIteration(options);
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

INSTANTIATE_TEST_SUITE_P(Generated, IteratorTest,
                        // proper indent
                        Combine(ValuesIn([]() -> std::set<autopas::ContainerOption> {
                                  auto allContainerOptions = autopas::allContainerOptions;
                                  /// @TODO no verletClusterLists yet, so we erase it for now.
                                  allContainerOptions.erase(
                                      allContainerOptions.find(autopas::ContainerOption::verletClusterLists));
                                  return allContainerOptions;
                                }()),
                                Values(0.5, 1., 1.5)),
                        IteratorTest::PrintToStringParamName());
