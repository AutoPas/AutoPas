/**
 * @file ParticleCounterTest.cpp
 * @author D. Martin
 * @date 10.09.23
 */

#include "ParticleCounterTest.h"

#include "autopas/tuning/selectors/ContainerSelector.h"

/**
 * Lambda to generate a readable string out of the parameter of this test.
 */
static auto toString = [](const auto &info) {
  auto containerOption = info.param;
  std::stringstream resStream;
  resStream << containerOption.to_string();
  std::string res = resStream.str();
  return res;
};

/**
 * Tests for all containers if the method getNumberOfParticles(IteratorBehavior behavior) returns the correct number of
 * particles internally stored in cells with respect to different IteratorBehaviors.
 *
 */
TEST_P(ParticleCounterTest, testGetNumberOfParticles) {
  const auto containerOption = GetParam();

  // Construct container
  constexpr double skin = _cutoff * 0.1;

  const auto containerInfo = autopas::ContainerSelectorInfo{_boxMin,
                                                            _boxMax,
                                                            _cutoff,
                                                            _cellSizeFactor,
                                                            skin,
                                                            32,
                                                            8,
                                                            autopas::LoadEstimatorOption::none,
                                                            _orderCellsByMortonIndex,
                                                            _useOptimizedLJFunctor,
                                                            _useCompactSoA,
                                                            _reserveVLSizes,
                                                            _bucketSortParticles,
                                                            _sortVerletLists,
                                                            _sortingFrequency};
  auto container = autopas::ContainerSelector<Molecule>::generateContainer(containerOption, containerInfo);

  // add owned particles
  for (int i = 0; i < _numOwnedParticles; i++) {
    Molecule p({1., 1., 1.}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    container->addParticle(p);
  }

  // add halo particles
  for (int i = 0; i < _numHaloParticles; i++) {
    Molecule p({_boxMax[0] + 1., _boxMax[1] + 1., _boxMax[2] + 1.}, {0., 0., 0.}, _numOwnedParticles + i);
    p.setOwnershipState(autopas::OwnershipState::halo);
    container->addHaloParticle(p);
  }

  // LinkedCellsReferences remove dummies during adding particles
  if (containerOption != autopas::options::ContainerOption::linkedCellsReferences) {
    // add dummy particle
    for (int i = 0; i < _numDummyParticles; i++) {
      Molecule p({1., 1., 1.}, {0., 0., 0.}, _numOwnedParticles + _numHaloParticles + i);
      p.setOwnershipState(autopas::OwnershipState::dummy);
      container->addParticle(p);
    }
  }

  // by default getNumberOfParticles() should return the number of owned particles
  EXPECT_EQ(container->getNumberOfParticles(), _numOwnedParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::owned), _numOwnedParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::halo), _numHaloParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::ownedOrHalo),
            _numOwnedParticles + _numHaloParticles);
  if (containerOption != autopas::options::ContainerOption::linkedCellsReferences) {
    EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::dummy), _numDummyParticles);
    EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::ownedOrHaloOrDummy),
              _numOwnedParticles + _numHaloParticles + _numDummyParticles);
  }
}

INSTANTIATE_TEST_SUITE_P(Generated, ParticleCounterTest, ::testing::ValuesIn(autopas::ContainerOption::getAllOptions()),
                         toString);