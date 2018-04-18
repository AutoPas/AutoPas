/**
 * @file Newton3OnOffTest.cpp
 * @author seckler
 * @date 18.04.18
 */

#include "Newton3OnOffTest.h"

using ::testing::_;       // anything is ok
using ::testing::Return;  // anything is ok

double Newton3OnOffTest::fRand(double fMin, double fMax) const {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

std::array<double, 3> Newton3OnOffTest::randomPosition(
    const std::array<double, 3> &boxMin,
    const std::array<double, 3> &boxMax) const {
  std::array<double, 3> r{};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}

void Newton3OnOffTest::fillContainerWithMolecules(
    unsigned long numMolecules,
    AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>>
        &autoPas) const {
  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(getBoxMin()), boxMax(getBoxMax());

  for (int i = 0; i < numMolecules; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::Particle m(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    autoPas.addParticle(m);
  }
}

TEST_F(Newton3OnOffTest, testAoS) {
  for (auto containerOption : autopas::possibleContainerOptions) {
    autoPas.init(getBoxMin(), getBoxMax(), getCutoff(), containerOption);
    fillContainerWithMolecules(100, autoPas);

    // with newton 3:
    int callsNewton3 = 0;
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(false));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true))
        .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsNewton3++; }));
    autoPas.iteratePairwise(&mockFunctor, autopas::DataLayoutOption::aos);

    // without newton 3:
    int callsNonNewton3 = 0;
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(false));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true))
        .Times(0);  // disables newton3 variant
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, false))
        .WillRepeatedly(
            testing::InvokeWithoutArgs([&]() { callsNonNewton3++; }));
    autoPas.iteratePairwise(&mockFunctor, autopas::DataLayoutOption::aos);

    EXPECT_EQ(callsNewton3 * 2,
              callsNonNewton3);  // should be called exactly two times
  }
}

TEST_F(Newton3OnOffTest, testSoA) {
  for (auto containerOption : autopas::possibleContainerOptions) {
    autoPas.init(getBoxMin(), getBoxMax(), getCutoff(), containerOption);
    /// @todo write test for SoA
    // autoPas.iteratePairwise(&mockFunctor, autopas::DataLayoutOption::soa);
  }
}