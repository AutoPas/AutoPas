/**
 * @file EvidenceCollectionTest.cpp
 * @author F. Gratl
 * @date 05.07.23
 */

#include "EvidenceCollectionTest.h"

#include "autopas/tuning/searchSpace/EvidenceCollection.h"

/**
 * Fill the collection with evidence from two phases and check what getOptimalConfiguration() and
 * getLatestOptimalConfiguration() return.
 */
TEST_F(EvidenceCollectionTest, testGetOptimalConfiguration) {
  autopas::EvidenceCollection collection;

  const autopas::Evidence evidenceLC_C01{0, 0, 100};
  collection.addEvidence(_configurationLC_C01, evidenceLC_C01);
  const autopas::Evidence evidenceLC_C08{5, 0, 50};
  collection.addEvidence(_configurationLC_C08, evidenceLC_C08);
  // one evidence from a different tuning phase
  const autopas::Evidence evidenceLC_Sliced{10, 1, 10};
  collection.addEvidence(_configurationLC_Sliced, evidenceLC_Sliced);
  collection.addEvidence(_configurationLC_C01, {15, 1, 100});

  {
    const auto &[optimalConf, optimalEvidence] = collection.getOptimalConfiguration(0);
    EXPECT_EQ(optimalConf, _configurationLC_C08);
    EXPECT_EQ(optimalEvidence, evidenceLC_C08);
  }

  {
    const auto &[optimalConf, optimalEvidence] = collection.getOptimalConfiguration(1);
    EXPECT_EQ(optimalConf, _configurationLC_Sliced);
    EXPECT_EQ(optimalEvidence, evidenceLC_Sliced);
  }

  {
    const auto &[optimalConf, optimalEvidence] = collection.getLatestOptimalConfiguration();
    EXPECT_EQ(optimalConf, _configurationLC_Sliced);
    EXPECT_EQ(optimalEvidence, evidenceLC_Sliced);
  }
}