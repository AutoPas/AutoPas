/**
 * @file Newton3OnOffTest.cpp
 * @author seckler
 * @date 18.04.18
 */

#include "Newton3OnOffTest.h"

using ::testing::_;       // anything is ok
using ::testing::Return;  // anything is ok

TEST_F(Newton3OnOffTest, testAoS) {
  for (auto containerOption : autopas::allContainerOptions) {
    // @todo remove this when cluster list can iterate with newton3 active
    if (containerOption == autopas::ContainerOption::verletClusterLists) {
      continue;
    }

    // needs two samples per container because we test with and without newton 3 and never take time measurements
    autoPas.setBoxMin(getBoxMin());
    autoPas.setBoxMax(getBoxMax());
    autoPas.setCutoff(getCutoff());
    autoPas.setVerletSkin(getVerletSkin());
    autoPas.setVerletRebuildFrequency(getVerletRebuildFrequency());
    autoPas.setAllowedContainers({containerOption});
    autoPas.setAllowedTraversals(autopas::allTraversalOptions);
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
    autoPas.setAllowedNewton3Options(autopas::allNewton3Options);
    autoPas.setTuningInterval(100);
    autoPas.setSelectorStrategy(autopas::SelectorStrategy::fastestAbs);
    autoPas.setNumSamples(1);

    autoPas.init();
    autopas::MoleculeLJ defaultParticle;
    RandomGenerator::fillWithParticles(*autoPas.getContainer(), defaultParticle, 100);
    RandomGenerator::fillWithHaloParticles(*autoPas.getContainer(), defaultParticle,
                                           autoPas.getContainer()->getCutoff(), 10);

    EXPECT_CALL(mockFunctor, isRelevantForTuning()).Times(testing::AtLeast(1)).WillRepeatedly(Return(true));

    // with newton 3:
    int callsNewton3 = 0;
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(false));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNewton3++;
    }));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true)).Times(testing::AtLeast(1));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, false)).Times(0);  // disables newton3 variant
    autoPas.iteratePairwise(&mockFunctor);

    // without newton 3:
    int callsNonNewton3 = 0;
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(false));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, false)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNonNewton3++;
    }));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, false)).Times(testing::AtLeast(1));
    EXPECT_CALL(mockFunctor, AoSFunctor(_, _, true)).Times(0);  // disables newton3 variant
    autoPas.iteratePairwise(&mockFunctor);

    EXPECT_EQ(callsNewton3 * 2, callsNonNewton3);  // should be called exactly two times

    if (::testing::Test::HasFailure()) {
      std::cerr << "Failures for containeroption: " << containerOption << std::endl;
    }
  }
}

TEST_F(Newton3OnOffTest, testSoA) {
  for (auto containerOption : autopas::allContainerOptions) {
    if (containerOption == autopas::ContainerOption::verletLists ||
        containerOption == autopas::ContainerOption::verletListsCells ||
        containerOption == autopas::ContainerOption::verletClusterLists) {
      continue;
    }
    autoPas.setBoxMax(getBoxMin());
    autoPas.setBoxMax(getBoxMax());
    autoPas.setCutoff(getCutoff());
    autoPas.setVerletSkin(getVerletSkin());
    autoPas.setVerletRebuildFrequency(getVerletRebuildFrequency());
    autoPas.setAllowedContainers({containerOption});
    autoPas.setAllowedTraversals(autopas::allTraversalOptions);
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::soa});
    autoPas.setAllowedNewton3Options(autopas::allNewton3Options);
    autoPas.setTuningInterval(100);
    autoPas.setSelectorStrategy(autopas::SelectorStrategy::fastestAbs);
    autoPas.setNumSamples(1);

    autoPas.init();
    autopas::MoleculeLJ defaultParticle;
    RandomGenerator::fillWithParticles(*autoPas.getContainer(), defaultParticle, 100);
    RandomGenerator::fillWithHaloParticles(*autoPas.getContainer(), defaultParticle,
                                           autoPas.getContainer()->getCutoff(), 10);

    EXPECT_CALL(mockFunctor, isRelevantForTuning()).Times(testing::AtLeast(1)).WillRepeatedly(Return(true));

    // loader and extractor will be called, we don't care how often.
    EXPECT_CALL(mockFunctor, SoALoader(_, _)).Times(testing::AtLeast(1));
    EXPECT_CALL(mockFunctor, SoAExtractor(_, _)).Times(testing::AtLeast(1));

    // with newton 3:
    int callsNewton3SC = 0;    // same cell
    int callsNewton3Pair = 0;  // pair of cells
    EXPECT_CALL(mockFunctor, allowsNewton3()).WillOnce(Return(true));
    EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillOnce(Return(false));

    // single cell
    EXPECT_CALL(mockFunctor, SoAFunctor(_, true)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNewton3SC++;
    }));
    EXPECT_CALL(mockFunctor, SoAFunctor(_, true)).Times(testing::AtLeast(1));

    // pair of cells
    EXPECT_CALL(mockFunctor, SoAFunctor(_, _, true)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNewton3Pair++;
    }));
    EXPECT_CALL(mockFunctor, SoAFunctor(_, _, true)).Times(testing::AtLeast(1));

    autoPas.iteratePairwise(&mockFunctor);

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
    EXPECT_CALL(mockFunctor, SoAFunctor(_, false)).Times(testing::AtLeast(1));

    // pair of cells
    EXPECT_CALL(mockFunctor, SoAFunctor(_, _, false)).WillRepeatedly(testing::InvokeWithoutArgs([&]() {
      callsNonNewton3Pair++;
    }));
    EXPECT_CALL(mockFunctor, SoAFunctor(_, _, false)).Times(testing::AtLeast(1));
    autoPas.iteratePairwise(&mockFunctor);

    EXPECT_EQ(callsNewton3SC, callsNonNewton3SC) << "for containeroption: " << containerOption;
    // should be called exactly two times
    EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair) << "for containeroption: " << containerOption;
    // should be called exactly two times
  }
}
