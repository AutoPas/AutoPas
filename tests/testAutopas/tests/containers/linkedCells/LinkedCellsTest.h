/**
 * @file LinkedCellsTest.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/particles/Particle.h"
#include "testingHelpers/commonTypedefs.h"

template <class TestingType>
class LinkedCellsTest : public AutoPasTestBase {
 public:
  using LinkedCellsType = std::tuple_element_t<0, TestingType>;
  using keepListsValid = std::tuple_element_t<1, TestingType>;
  LinkedCellsTest() : _linkedCells({0., 0., 0.}, {10., 10., 10.}, 1., 0., 1.) {}

 protected:
  LinkedCellsType _linkedCells;
  bool _keepListsValid{keepListsValid()};

  void checkParticleIDsInCells(
      LinkedCellsType &linkedCells,
      std::map<unsigned long /*cellID*/, std::vector<int> /*particleIDs*/> expectedParticleIDsInCells, bool ordered,
      int line);
};
template <class TestingType>
void LinkedCellsTest<TestingType>::checkParticleIDsInCells(
    LinkedCellsType &linkedCells, std::map<unsigned long, std::vector<int>> expectedParticleIDsInCells, bool ordered,
    int line) {
  // check particles are where we expect them to be (and nothing else)
  for (size_t i = 0; i < linkedCells.getCells().size(); ++i) {
    if (auto iter = expectedParticleIDsInCells.find(i); iter != expectedParticleIDsInCells.end()) {
      auto expectedIDs = iter->second;
      auto expectedNumParticles = expectedIDs.size();
      ASSERT_EQ(linkedCells.getCells()[i].numParticles(), expectedNumParticles) << "called from line: " << line;

      std::vector<int> foundIDs;
      {
        auto pIter = linkedCells.getCells()[i].begin();
        for (int j = 0; j < expectedNumParticles; ++j, ++pIter) {
          ASSERT_TRUE(pIter.isValid()) << "called from line: " << line;
          foundIDs.push_back(pIter->getID());
        }
      }
      if (ordered) {
        EXPECT_THAT(foundIDs, testing::ElementsAreArray(expectedIDs)) << "called from line: " << line;
      } else {
        EXPECT_THAT(foundIDs, testing::UnorderedElementsAreArray(expectedIDs)) << "called from line: " << line;
      }
    } else {
      EXPECT_TRUE(linkedCells.getCells()[i].isEmpty()) << "Cell: " << i << std::endl << "called from line: " << line;
    }
  }
}
