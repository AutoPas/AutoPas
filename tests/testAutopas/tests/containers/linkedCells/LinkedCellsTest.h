
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
  /// @todo c++20 replace with:
  /// using LinkedCellsType = std::tuple_element_t<0, TestingType>;
  /// using keepListsValid = std::tuple_element_t<1, TestingType>;
  using LinkedCellsType = typename TestingType::first_t;
  using KeepListsValid = typename TestingType::second_t;
  LinkedCellsTest() {}

 protected:
  bool _keepListsValid{KeepListsValid()};

  void checkParticleIDsInCells(
      LinkedCellsType &linkedCells,
      std::map<unsigned long /*cellID*/, std::vector<std::tuple<int, autopas::OwnershipState>> /*particleIDs*/>
          expectedParticleIDsInCells,
      bool ordered, int line);
};

template <class TestingType>
void LinkedCellsTest<TestingType>::checkParticleIDsInCells(
    LinkedCellsType &linkedCells,
    std::map<unsigned long, std::vector<std::tuple<int, autopas::OwnershipState>>> expectedParticleIDsInCells,
    bool ordered, int line) {
  // check particles are where we expect them to be (and nothing else)
  for (size_t i = 0; i < linkedCells.getCells().size(); ++i) {
    if (auto iter = expectedParticleIDsInCells.find(i); iter != expectedParticleIDsInCells.end()) {
      auto expectedParticles = iter->second;
      auto expectedNumParticles = expectedParticles.size();
      ASSERT_EQ(linkedCells.getCells()[i].size(), expectedNumParticles) << "called from line: " << line;

      std::vector<std::tuple<int, autopas::OwnershipState>> foundParticles;
      {
        auto pIter = linkedCells.getCells()[i].begin();
        for (int j = 0; j < expectedNumParticles; ++j, ++pIter) {
          foundParticles.push_back({pIter->getID(), pIter->getOwnershipState()});
        }
      }
      if (ordered) {
        EXPECT_THAT(foundParticles, testing::ElementsAreArray(expectedParticles)) << "called from line: " << line;
      } else {
        EXPECT_THAT(foundParticles, testing::UnorderedElementsAreArray(expectedParticles))
            << "called from line: " << line;
      }
    } else {
      bool isEmpty = linkedCells.getCells()[i].isEmpty();
      if (not isEmpty) {
        // If the cell isn't empty, all particles should be dummy particles!
        for (auto &p : linkedCells.getCells()[i]) {
          EXPECT_TRUE(p.isDummy()) << "Cell: " << i << " Particle ID: " << p.getID() << std::endl
                                   << "called from line: " << line;
        }
      }
    }
  }
}
