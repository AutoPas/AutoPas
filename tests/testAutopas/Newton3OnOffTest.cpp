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

std::array<double, 3> Newton3OnOffTest::randomPosition(const std::array<double, 3> &boxMin,
                                                       const std::array<double, 3> &boxMax) const {
  std::array<double, 3> r{};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}

void Newton3OnOffTest::fillContainerWithMolecules(
    unsigned long numMolecules,
    AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> &autoPas) const {
  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(getBoxMin()), boxMax(getBoxMax());

  for (size_t i = 0; i < numMolecules; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::Particle m(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
    autoPas.addParticle(m);
  }
}

TEST_F(Newton3OnOffTest, testAoS) {
  for (auto containerOption : autopas::allContainerOptions) {
    autoPas.init(getBoxMin(), getBoxMax(), getCutoff(), {containerOption});
    fillContainerWithMolecules(100, autoPas);

    // with newton 3:
    int callsNewton3 = 0;
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(false));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNewton3++;
    }));
    autoPas.iteratePairwise(&mockFunctor, autopas::DataLayoutOption::aos);

    // without newton 3:
    int callsNonNewton3 = 0;
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(false));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true)).Times(0);  // disables newton3 variant
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, false)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNonNewton3++;
    }));
    autoPas.iteratePairwise(&mockFunctor, autopas::DataLayoutOption::aos);

    EXPECT_EQ(callsNewton3 * 2,
              callsNonNewton3);  // should be called exactly two times
  }
}

TEST_F(Newton3OnOffTest, testSoA) {
  for (auto containerOption : autopas::allContainerOptions) {
    autoPas.init(getBoxMin(), getBoxMax(), getCutoff(), {containerOption});
    fillContainerWithMolecules(100, autoPas);

    // loader and extractor will be called, we don't care how often.
    EXPECT_CALL(mockFunctor, SoALoader(_, _, _)).Times(testing::AtLeast(1));
    EXPECT_CALL(mockFunctor, SoAExtractor(_, _, _)).Times(testing::AtLeast(1));

    // with newton 3:
    int callsNewton3SC = 0;    // same cell
    int callsNewton3Pair = 0;  // pair of cells
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(false));

    // single cell
    EXPECT_CALL(mockFunctor, SoAFunctor(_, true)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNewton3SC++;
    }));

    // pair of cells
    EXPECT_CALL(mockFunctor, SoAFunctor(_, _, true)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNewton3Pair++;
    }));

    autoPas.iteratePairwise(&mockFunctor, autopas::DataLayoutOption::soa);

    // without newton 3:
    int callsNonNewton3SC = 0;
    int callsNonNewton3Pair = 0;
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(false));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, SoAFunctor(_, _, true)).Times(0);  // disables newton3 variant

    // single cell
    EXPECT_CALL(mockFunctor, SoAFunctor(_, false)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNonNewton3SC++;
    }));

    // pair of cells
    EXPECT_CALL(mockFunctor, SoAFunctor(_, _, false)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNonNewton3Pair++;
    }));
    autoPas.iteratePairwise(&mockFunctor, autopas::DataLayoutOption::soa);

    EXPECT_EQ(callsNewton3SC,
              callsNonNewton3SC);  // should be called exactly two times
    EXPECT_EQ(callsNewton3Pair * 2,
              callsNonNewton3Pair);  // should be called exactly two times
  }
}