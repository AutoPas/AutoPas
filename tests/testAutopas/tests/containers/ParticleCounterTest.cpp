/**
 * @file ParticleCounterTest.cpp
 * @author D. Martin
 * @date 10.09.23
 */

#include "ParticleCounterTest.h"

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

TEST_P(ParticleCounterTest, testGetNumberOfParticles) {
  auto containerOption = GetParam();

  // Construct container
  autopas::ContainerSelector<Molecule> selector{_boxMin, _boxMax, _cutoff};
  constexpr double skinPerTimestep = _cutoff * 0.1;
  constexpr unsigned int rebuildFrequency = 1;
  selector.selectContainer(containerOption,
                           autopas::ContainerSelectorInfo{_cellSizeFactor, skinPerTimestep, rebuildFrequency, 32,
                                                          autopas::LoadEstimatorOption::none});
  auto &container = selector.getCurrentContainer();

  // add owned particles
  for (int i = 0; i < _numOwnedParticles; i++) {
    Molecule p({1., 1., 1.}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    container.addParticle(p);
  }

  // add halo particles
  for (int i = 0; i < _numHaloParticles; i++) {
    Molecule p({_boxMax[0] + 1., _boxMax[1] + 1., _boxMax[2] + 1.}, {0., 0., 0.}, _numOwnedParticles + i);
    p.setOwnershipState(autopas::OwnershipState::halo);
    container.addHaloParticle(p);
  }

  // LinkedCellsReferences remove dummies during adding particles
  if (not(containerOption == autopas::options::ContainerOption::linkedCellsReferences)) {
    // add dummy particle
    for (int i = 0; i < _numDummyParticles; i++) {
      Molecule p({1., 1., 1.}, {0., 0., 0.}, _numOwnedParticles + _numHaloParticles + i);
      p.setOwnershipState(autopas::OwnershipState::dummy);
      container.addParticle(p);
    }
  }

  // by default getNumberOfParticles() should return the number of owned particles
  EXPECT_EQ(container.getNumberOfParticles(), _numOwnedParticles);
  EXPECT_EQ(container.getNumberOfParticles(autopas::options::IteratorBehavior::owned), _numOwnedParticles);
  EXPECT_EQ(container.getNumberOfParticles(autopas::options::IteratorBehavior::halo), _numHaloParticles);
  EXPECT_EQ(container.getNumberOfParticles(autopas::options::IteratorBehavior::ownedOrHalo),
            _numOwnedParticles + _numHaloParticles);
  if (not(containerOption == autopas::options::ContainerOption::linkedCellsReferences)) {
    EXPECT_EQ(container.getNumberOfParticles(autopas::options::IteratorBehavior::dummy), _numDummyParticles);
    EXPECT_EQ(container.getNumberOfParticles(autopas::options::IteratorBehavior::ownedOrHaloOrDummy),
              _numOwnedParticles + _numHaloParticles + _numDummyParticles);
  }
}

INSTANTIATE_TEST_SUITE_P(Generated, ParticleCounterTest, ::testing::ValuesIn(autopas::ContainerOption::getAllOptions()),
                         toString);