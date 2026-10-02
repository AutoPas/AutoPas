/**
 * @file EvidenceCollectionTest.cpp
 * @author F. Gratl
 * @date 05.07.23
 */

#include "EvidenceCollectionTest.h"

#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "testingHelpers/ArbitraryConfigurations.h"

/**
 * Fill the collection with evidence from two phases and check what getOptimalConfiguration() and
 * getLatestOptimalConfiguration() return.
 */
TEST_F(EvidenceCollectionTest, testGetOptimalConfiguration) {
  autopas::EvidenceCollection collection;

  const autopas::Evidence evidence0{0, 0, 100};
  collection.addEvidence(arbitraryConfigurations::_arbitrary_config_2B_0, evidence0);
  const autopas::Evidence evidence1{5, 0, 50};
  collection.addEvidence(arbitraryConfigurations::_arbitrary_config_2B_1, evidence1);
  // one evidence from a different tuning phase
  const autopas::Evidence evidence2{10, 1, 10};
  collection.addEvidence(arbitraryConfigurations::_arbitrary_config_2B_2, evidence2);
  collection.addEvidence(arbitraryConfigurations::_arbitrary_config_2B_3, {15, 1, 100});

  {
    const auto &[optimalConf, optimalEvidence] = collection.getOptimalConfiguration(0);
    EXPECT_EQ(optimalConf, arbitraryConfigurations::_arbitrary_config_2B_1);
    EXPECT_EQ(optimalEvidence, evidence1);
  }

  {
    const auto &[optimalConf, optimalEvidence] = collection.getOptimalConfiguration(1);
    EXPECT_EQ(optimalConf, arbitraryConfigurations::_arbitrary_config_2B_2);
    EXPECT_EQ(optimalEvidence, evidence2);
  }

  {
    const auto &[optimalConf, optimalEvidence] = collection.getLatestOptimalConfiguration();
    EXPECT_EQ(optimalConf, arbitraryConfigurations::_arbitrary_config_2B_2);
    EXPECT_EQ(optimalEvidence, evidence2);
  }
}